#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>
#include <gmx_reader.h>
#include "HBondDistribution.h"

#define PI 3.14159265359 

using namespace std;

// First define some new class member functions
// Default constructor
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
model::model( string _inpf_ ) : gmx_reader::gmx_reader( _inpf_ )
{
// set variables for this subclass
    for ( int i = 0; i < nuParams; i ++ )
    {
        if ( uParams[i] == "db" )       db      = stof( uValues[i] );
        if ( uParams[i] == "dr" )       dr      = stof( uValues[i] );
        if ( uParams[i] == "rmin" )     r_min   = stof( uValues[i] );
        if ( uParams[i] == "rmax" )     r_max   = stof( uValues[i] );
        if ( uParams[i] == "bmin" )     b_min   = stof( uValues[i] );
        if ( uParams[i] == "bmax" )     b_max   = stof( uValues[i] );
        if ( uParams[i] == "gbrfname" )  gbrfname = uValues[i];
        if ( uParams[i] == "gbrthbfname" )  gbrthbfname = uValues[i];
        if ( uParams[i] == "pmffname" )  pmffname = uValues[i];
        if ( uParams[i] == "tcffname" )  tcffname = uValues[i];
        if ( uParams[i] == "rhbond" )   rhbond   = stof(uValues[i]);
        if ( uParams[i] == "betahbond") betahbond= stof(uValues[i]);
        if ( uParams[i] == "deltaTCFMax" ) deltaTCFMax = stoi(uValues[i]);
    }

    cout << "Set db to: "   << db << endl;
    cout << "Set dr to: "   << dr << endl;
    cout << "Set rmin to: " << r_min << endl;
    cout << "Set rmax to: " << r_max << endl;
    cout << "Set bmin to: " << b_min << endl;
    cout << "Set bmax to: " << b_max << endl;
    cout << "Set gbrfname to: " << gbrfname << endl;
    cout << "Set gbrthbfname to: " << gbrthbfname << endl;
    cout << "Set pmffname to: " << pmffname << endl;
    cout << "Set tcffname to: " << tcffname << endl;
    cout << "Set rhbond to: " << rhbond << endl;
    cout << "Set betahbond to: " << betahbond << endl;
    cout << "Set deltaTCFMax to: " << deltaTCFMax << endl;

    //  make sure deltaTCFMax is okay -- adjust if you need to
    deltaTCFMax = min(deltaTCFMax,nsamples);

    // allocate space for nearest neighbor classifications
    npoints_r    = (int) round( (r_max-r_min)/dr + 1);           // total number of grid points
    npoints_b    = (int) round( (b_max-b_min)/db + 1);           // total number of grid points

    // the array for the H-bond distribution function
    gbr         = new double[npoints_b*npoints_r]();
    pmf         = new double[npoints_b*npoints_r]();
    hbonded     = new bool[nmol*nmol]();
    hbonded_t0  = new bool[nmol*nmol]();
    hbonded_t   = new bool[nmol*nmol]();
    hbondTCF    = new float[deltaTCFMax]();
    rr          = new float[ nmol*nmol ]();
    bb1         = new float[ nmol*nmol ]();
    bb2         = new float[ nmol*nmol ]();
    gbr_thb     = new double[npoints_b*npoints_r]();

}

// Default Destructor
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
model::~model()
{
    delete [] gbr;
    delete [] pmf;
    delete [] hbonded;
    delete [] hbonded_t0;
    delete [] hbonded_t;
    delete [] hbondTCF;
    delete [] rr;
    delete [] bb1;
    delete [] bb2;
    delete [] gbr_thb;
}

float model::get_r( int mol1, int mol2 )
{
    // Return the OO distance between two molecules
    float OOvec[3];
    int i;

    // Get the OO vector
    for ( i = 0; i < 3; i ++ ) OOvec[i] = x[ mol2 * natoms_mol ][i] \
                                        - x[ mol1 * natoms_mol ][i];
    // return the OO vector magnitude
    minImage( OOvec );
    return mag3( OOvec );
}

float model::get_b( int mol1, int mol2, int whichH )
{
    float OOvec[3], OHvec[3];
    int i;
    float beta;
    float arg;

    if ( whichH < 1 or whichH > 2 ){
        cout << "Bad value for whichH: " << whichH << endl;
        exit(1);
    }

    // Return the OH-O angle between molecule 1 and 2
    // considers molecule 1 only as the H-bond acceptor
    for ( i = 0; i < 3; i ++ ) OOvec[i] = x[ mol2 * natoms_mol ][i] \
                                        - x[ mol1 * natoms_mol ][i];
    minImage( OOvec );

    // get the OH vectors for the reference molecule
    for ( i = 0; i < 3; i ++ ) OHvec[i] = x[ mol1 * natoms_mol + whichH ][i] \
                                        - x[ mol1 * natoms_mol ][i];
    minImage( OHvec );

    // determine beta
    arg   = dot3(OOvec,OHvec)/(mag3(OOvec)*mag3(OHvec));
    // note acos is only defined when the argument is between zero and 1 -- this next line makes sure
    // the program doesnt break due to numerical errors that cause the argument to go slightly over 1 or slightly below -1
    arg   = max( arg, (float) -1. ); arg = min( arg, (float) 1. );
    if ( arg < -1 or arg > 1. ) cout << "WARNING: arg out of bounds" << endl;
    beta  = 180./PI*acos(arg);

    // return beta
    return  beta;
}

int model::get_rnx( float r )
{
    return round((r-r_min)/dr);
}

int model::get_bnx( float b )
{
    return round((b-b_min)/db);
}
                
int model::get_nx( int rnx, int bnx )
{
    return bnx*npoints_r + rnx;
}

void model::write_gbr()
{
    int rnx, bnx, nx;

    FILE *file = fopen(gbrfname.c_str(),"w");
    fprintf( file, "#R (nm) Beta (deg) g(R,Beta)\n");

    for ( rnx = 0; rnx < npoints_r; rnx ++ ){
        for ( bnx = 0; bnx < npoints_b; bnx ++ ){
            nx = get_nx( rnx, bnx );
            fprintf( file, "%g %g %g\n", r_min + rnx*dr + dr/2., b_min + bnx*db + db/2., gbr[ nx ]);
        }
    }
    fclose( file );
}

void model::write_pmf()
{
    int rnx, bnx, nx;

    FILE *file = fopen(pmffname.c_str(),"w");
    fprintf( file, "#R (nm) Beta (deg) PMF(R,Beta)\n");

    for ( rnx = 0; rnx < npoints_r; rnx ++ ){
        for ( bnx = 0; bnx < npoints_b; bnx ++ ){
            nx = get_nx( rnx, bnx );
            fprintf( file, "%g %g %g\n", r_min + rnx*dr + dr/2., b_min + bnx*db + db/2., pmf[ nx ]);
        }
    }
    fclose( file );
}

int model::getarraynx( int mol1, int mol2 )
{
    return mol1 * nmol + mol2;
}

bool model::is_hbond( float r, float beta1, float beta2 )
{
    float beta;

    // hbond condition is r<rhbond, beta<betahbond
    if ( r > rhbond ) return false;

    // only the minimum beta will be relevant here
    beta = min(beta1, beta2);
    if ( beta > betahbond ) return false;

    return true;
}

void model::write_hbond(int currentSample)
{
    string fname = ".hbond-" + to_string(currentSample) + ".dat";
    FILE *file = fopen( fname.c_str(),"wb");
    fwrite( hbonded, sizeof(bool), nmol*nmol, file );
    fwrite( rr , sizeof(float), nmol*nmol, file );
    fwrite( bb1, sizeof(float), nmol*nmol, file );
    fwrite( bb2, sizeof(float), nmol*nmol, file );
    fclose(file);
}

void model::read_hbond_t0(int currentSample)
{
    string fname = ".hbond-" + to_string(currentSample) + ".dat";
    FILE *file = fopen( fname.c_str(),"rb");
    fread( hbonded_t0, sizeof(bool), nmol*nmol, file );
    fread( rr  , sizeof(float), nmol*nmol, file );
    fread( bb1 , sizeof(float), nmol*nmol, file );
    fread( bb2 , sizeof(float), nmol*nmol, file );
    fclose( file );
}

void model::read_hbond_t(int currentSample)
{
    string fname = ".hbond-" + to_string(currentSample) + ".dat";
    FILE *file = fopen( fname.c_str(),"rb");
    fread( hbonded_t, sizeof(bool), nmol*nmol, file );
    fread( rr  , sizeof(float), nmol*nmol, file );
    fread( bb1 , sizeof(float), nmol*nmol, file );
    fread( bb2 , sizeof(float), nmol*nmol, file );
    fclose( file );
}

void model::remove_hbond_files()
{
    int currentSample;
    string fname;

    for ( currentSample = 0; currentSample < nsamples; currentSample ++ ){
        fname = ".hbond-" + to_string(currentSample) + ".dat";
        remove(fname.c_str());
    }

}

float model::get_hbond_TCF( int deltaTCF )
{
    float tcf = 0;
    int sample, nTCFsamples;
    int i, mol1, mol2;
    int rnx, bnx, nx, arraynx;
    int kount;
    double ht = 0.;
    float bb;

    kount = 0;
    nTCFsamples = nsamples - deltaTCF; 

    for ( i = 0; i < npoints_r*npoints_b; i ++ ) gbr_thb[i] = 0.;

    for ( sample = 0; sample < nTCFsamples; sample ++ )
    {
        read_hbond_t0(sample);
        read_hbond_t(sample+deltaTCF);

        // hydrogen bond time correlation function
        for ( i = 0; i < nmol*nmol; i ++ ){
            ht += hbonded_t0[i]*hbonded_t[i];
        }

        // gbr(t|hb)
        for ( mol1 = 0; mol1 < nmol; mol1 ++ ){
            for ( mol2 = 0; mol2 < nmol; mol2 ++ ){
                
                if ( mol1 == mol2 ) continue;

                kount += 2;

                // get index for mol1 and mol2
                arraynx = getarraynx( mol1, mol2 );

                // get r
                if ( rr [ arraynx ] > r_max or rr[ arraynx ] < r_min ) continue;
                rnx = get_rnx( rr[ arraynx ] );

                // it seems like in the other paper they consider only the beta of the hydrogen that is initially hydrogen bonded
                // WARNING -- THIS MAY NOT QUITE BE RIGHT -- YOU MAY NEED TO THINK ABOUT THIS A LITTLE MORE
                // It may be that you need to define hbonded_t0 not in terms of molecules, but in terms of specific OHO pairs.
                bb = min( bb1[ arraynx ], bb2[ arraynx ] );
                if ( bb <= b_max and bb >= b_min ){

                    bnx = get_bnx( bb );
                    nx  = get_nx( rnx, bnx );
                    gbr_thb[ nx ] += 1. * hbonded_t0[ arraynx ];
                }
            }
        }
    }

    // Normalization
    ht /= 1.*nTCFsamples;
    for ( i = 0; i < npoints_r*npoints_b; i ++ ) gbr_thb[i] /= 1.*kount*gbr[i];

    // may include option to write only certian ones that you want, but for now this is ok
    write_gbr_thb( deltaTCF );

    return ht;
}

void model::write_hbond_tcf()
{
    int deltaTCF;

    FILE *file = fopen(tcffname.c_str(),"w");
    fprintf( file, "#t (ns) hbond TCF\n");

    for ( deltaTCF = 0; deltaTCF < deltaTCFMax; deltaTCF ++ )
    {
        fprintf( file, "%g %g \n", deltaTCF*sampleEvery, hbondTCF[deltaTCF] );
    }
    fclose( file );
}

void model::write_gbr_thb( int deltaTCF )
{
    int rnx, bnx, nx;

    string fname = gbrthbfname+"-"+to_string(deltaTCF)+".dat";
    FILE *file = fopen(fname.c_str(),"w");
    fprintf( file, "#Time: %g\n", deltaTCF*sampleEvery );
    fprintf( file, "#R (nm) Beta (deg) g(R,Beta,t)\n");

    for ( rnx = 0; rnx < npoints_r; rnx ++ ){
        for ( bnx = 0; bnx < npoints_b; bnx ++ ){
            nx = get_nx( rnx, bnx );
            fprintf( file, "%g %g %g\n", r_min + rnx*dr + dr/2., b_min + bnx*db + db/2., gbr_thb[ nx ]);
        }
    }
    fclose( file );
}



// ************************************************************************ 
int main( int argc, char* argv[] )
{

    int currentSample, i, mol1, mol2, rnx, bnx, nx, arraynx;
    long int kount;
    float r, b1, b2;
    float ht0;

    // Check program input
    if ( argc != 2 ){
        printf("Program expects only one argument, which is the name of \n\
                an input file containing the details of the analysis.\nAborting...\n");
        exit(EXIT_FAILURE);
    }

    // get filename
    string inpf(argv[1]);

    // attempt to initialize reading of gmx file
    model reader( inpf );

    // count the number of points taken
    kount = 0;

    // loop over trajectory
    for ( currentSample = 0; currentSample < reader.nsamples; currentSample ++ ){
        cout << "\rCurrent time: " << reader.gmxtime << setprecision(2) << fixed <<  " (ps)";
        cout.flush();

        // reset hbond array
        memset( reader.hbonded, false, reader.nmol*reader.nmol*sizeof(bool));
        
        // calculate H-bond distribution of angles and distances 
        // loop over all pairs of molecules
        for ( mol1 = 0; mol1 < reader.nmol; mol1 ++ ){
            for ( mol2 = 0; mol2 < reader.nmol; mol2 ++ ){

                if ( mol1 == mol2 ) continue;

                // increment counter -- must divide at the end to get a probability -- each water has 2 OH angles
                kount += 2;

                // Get OO distance
                arraynx = reader.getarraynx( mol1, mol2 );
                // mol1 * reader.nmol + mol2
                reader.rr[ arraynx ] = reader.get_r( mol1, mol2 );

                // Get Angle -- Hydrogen 1
                reader.bb1[ arraynx ] = reader.get_b( mol1, mol2, 1 );
                if ( reader.bb1[ arraynx ] < reader.b_max and reader.bb1[ arraynx ] > reader.b_min ){
                    if ( reader.rr[ arraynx ] <= reader.r_max and reader.rr[ arraynx ] >= reader.r_min ){
                        rnx = reader.get_rnx( reader.rr[ arraynx ] );
                        bnx = reader.get_bnx( reader.bb1[arraynx ] );
                        nx  = reader.get_nx( rnx, bnx );
                        reader.gbr[ nx ] += 1.;
                    }
                }

                // Get Angle -- Hydrogen 2
                reader.bb2[ arraynx ] = reader.get_b( mol1, mol2, 2 );
                if ( reader.bb2[ arraynx ] < reader.b_max and reader.bb2[ arraynx ] > reader.b_min ){
                    if ( reader.rr[ arraynx ] <= reader.r_max and reader.rr[ arraynx ] >= reader.r_min ){
                        rnx = reader.get_rnx( reader.rr[ arraynx ] );
                        bnx = reader.get_bnx( reader.bb2[arraynx ] );
                        nx  = reader.get_nx( rnx, bnx );
                        reader.gbr[ nx ] += 1.;
                    }
                }

                // determine if mol1 and mol2 are hbonded
                if ( reader.is_hbond( reader.rr[ arraynx ], reader.bb1[ arraynx ], reader.bb2[ arraynx] ) ){
                    reader.hbonded[ arraynx ] = true;
                }
            }
        }
        // write h-bond matrix, r matrix, and beta matrix to a scratch file
        reader.write_hbond(currentSample);

        // Advance to next frame if we need that frame
        if ( currentSample != reader.nsamples - 1 ) reader.search_for_sample( currentSample + 1 );
    } 

    // normalize the probability distribution function and calculate the PMF
    for ( rnx = 0; rnx < reader.npoints_r; rnx ++ ){
        for ( bnx = 0; bnx < reader.npoints_b; bnx ++ ){
            nx = reader.get_nx( rnx, bnx );
            reader.gbr[ nx ] /= kount*1.;

            // note that the PMF is normalized by kT here
            reader.pmf[ nx ] = -1.*log(reader.gbr[nx]);
        }
    }

    // now calculate the h-bond correlation function
    for ( int deltaSample = 0; deltaSample < reader.deltaTCFMax; deltaSample ++ )
    {
        cout << "\rNow calculating hbond TCF at t= " << deltaSample*reader.sampleEvery << " (ps)";
        cout.flush();
        reader.hbondTCF[deltaSample] = reader.get_hbond_TCF( deltaSample );
        if ( deltaSample == 0 ){
            ht0 = reader.hbondTCF[0];
        }
        reader.hbondTCF[deltaSample] /= ht0;

    }
    // remove hbond matrix scratch files
    reader.remove_hbond_files();

    // write the probability distribution function and pmfs to a file
    reader.write_gbr();
    reader.write_pmf();
    reader.write_hbond_tcf();

    cout << endl << "DONE!" << endl;
}
