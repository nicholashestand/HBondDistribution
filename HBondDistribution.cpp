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
        if ( uParams[i] == "db" )               db              = stof( uValues[i] );
        if ( uParams[i] == "dr" )               dr              = stof( uValues[i] );
        if ( uParams[i] == "rmin" )             r_min           = stof( uValues[i] );
        if ( uParams[i] == "rmax" )             r_max           = stof( uValues[i] );
        if ( uParams[i] == "bmin" )             b_min           = stof( uValues[i] );
        if ( uParams[i] == "bmax" )             b_max           = stof( uValues[i] );
        if ( uParams[i] == "gbrfname" )         gbrfname        = uValues[i];
        if ( uParams[i] == "gbrthbfname" )      gbrthbfname     = uValues[i];
        if ( uParams[i] == "pmffname" )         pmffname        = uValues[i];
        if ( uParams[i] == "tcffname" )         tcffname        = uValues[i];
        if ( uParams[i] == "phbfname" )         phbfname        = uValues[i];
        if ( uParams[i] == "prfname" )          prfname         = uValues[i];
        if ( uParams[i] == "ptfname" )          ptfname         = uValues[i];
        if ( uParams[i] == "rhbond_max" )       rhbond_max      = stof(uValues[i]);
        if ( uParams[i] == "rhbond_min" )       rhbond_min      = stof(uValues[i]);
        if ( uParams[i] == "rhbond_maxT" )      rhbond_maxT     = stof(uValues[i]);
        if ( uParams[i] == "betahbond_max")     betahbond_max   = stof(uValues[i]);
        if ( uParams[i] == "betahbond_min")     betahbond_min   = stof(uValues[i]);
        if ( uParams[i] == "betahbond_maxR")    betahbond_maxR  = stof(uValues[i]);
        if ( uParams[i] == "deltaTCFMax" )      deltaTCFMax     = stoi(uValues[i]);
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
    cout << "Set rhbond_max to: " << rhbond_max << endl;
    cout << "Set rhbond_min to: " << rhbond_min << endl;
    cout << "Set rhbond_maxT to: " << rhbond_maxT << endl;
    cout << "Set betahbond_max to: " << betahbond_max << endl;
    cout << "Set betahbond_min to: " << betahbond_min << endl;
    cout << "Set betahbond_maxR to: " << betahbond_maxR << endl;
    cout << "Set deltaTCFMax to: " << deltaTCFMax << endl;

    //  make sure deltaTCFMax is okay -- adjust if you need to
    deltaTCFMax = min(deltaTCFMax,nsamples);

    // allocate space for nearest neighbor classifications
    npoints_r    = (int) round( (r_max-r_min)/dr + 1);           // total number of grid points
    npoints_b    = (int) round( (b_max-b_min)/db + 1);           // total number of grid points

    // the array for the H-bond distribution function
    gbr         = new double[npoints_b*npoints_r]();
    pmf         = new double[npoints_b*npoints_r]();
    hbonded     = new bool[nmol*nmol*2]();
    hbonded_t0  = new bool[nmol*nmol*2]();
    hbonded_t   = new bool[nmol*nmol*2]();
    hbondTCF    = new double[deltaTCFMax]();
    NHBt        = new double[deltaTCFMax]();
    NRt         = new double[deltaTCFMax]();
    NTt         = new double[deltaTCFMax]();
    rr          = new float[ nmol*nmol ]();
    bb          = new float[ nmol*nmol*2 ]();
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
    delete [] bb;
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

int model::getrrnx( int mol1, int mol2 )
{
    return mol1 * nmol + mol2;
}

int model::getbbnx( int mol1, int mol2, int h )
{
    return (mol1 * nmol + mol2)*2 + h-1;
}

bool model::is_hbond( float r, float beta )
{
    // hbond condition is r<rhbond, beta<betahbond
    if ( r >= rhbond_max or r <= rhbond_min ) return false;
    if ( beta >= betahbond_max or beta <= betahbond_min ) return false;
    return true;
}

void model::write_hbond(int currentSample)
{
    string fname = ".hbond-" + to_string(currentSample) + ".dat";
    FILE *file = fopen( fname.c_str(),"wb");
    fwrite( hbonded, sizeof(bool), nmol*nmol*2, file );
    fwrite( rr , sizeof(float), nmol*nmol, file );
    fwrite( bb, sizeof(float), nmol*nmol*2, file );
    fclose(file);
}

void model::read_hbond_t0(int currentSample)
{
    string fname = ".hbond-" + to_string(currentSample) + ".dat";
    FILE *file = fopen( fname.c_str(),"rb");
    fread( hbonded_t0, sizeof(bool), nmol*nmol*2, file );
    fread( rr  , sizeof(float), nmol*nmol, file );
    fread( bb , sizeof(float), nmol*nmol*2, file );
    fclose( file );
}

void model::read_hbond_t(int currentSample)
{
    string fname = ".hbond-" + to_string(currentSample) + ".dat";
    FILE *file = fopen( fname.c_str(),"rb");
    fread( hbonded_t, sizeof(bool), nmol*nmol*2, file );
    fread( rr  , sizeof(float), nmol*nmol, file );
    fread( bb , sizeof(float), nmol*nmol*2, file );
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

void model::get_hbond_dynamics( int deltaTCF )
{
    int sample, nTCFsamples, i;
    int rnx, bnx, nx, rrnx, bbnx, mol1, mol2;
    float beta, r;

    // Determine number of samples to take
    nTCFsamples = nsamples - deltaTCF; 

    // Initialize variables
    hbondTCF[ deltaTCF ] = 0.;
    NHBt[ deltaTCF ] = 0.;
    NRt[ deltaTCF ] = 0.;
    NRt[ deltaTCF ] = 0.;
    for ( i = 0; i < npoints_r*npoints_b; i ++ ) gbr_thb[i] = 0.;

    for ( sample = 0; sample < nTCFsamples; sample ++ )
    {
        // read hbond info from files
        read_hbond_t0(sample);
        read_hbond_t(sample+deltaTCF);

        // hydrogen bond time correlation function
        // NOTE: this considers hbonds between O-HO pairs -- not between pairs of molecules
        for ( i = 0; i < nmol*nmol*2; i ++ ){
            hbondTCF[ deltaTCF ] += hbonded_t0[i]*hbonded_t[i];
        }

        // gbr(t|hb)
        for ( mol1 = 0; mol1 < nmol; mol1 ++ ){
            for ( mol2 = 0; mol2 < nmol; mol2 ++ ){
                
                if ( mol1 == mol2 ) continue;

                // get rr index for mol1 and mol2
                rrnx = getrrnx( mol1, mol2 );

                // get r
                r = rr[ rrnx ];
                if ( r > r_max or r < r_min ) continue;
                rnx = get_rnx( r );

                // hydrogen 1
                bbnx = getbbnx( mol1, mol2, 1 );
                beta = bb[ bbnx ];
                bnx = get_bnx( beta );
                nx  = get_nx( rnx, bnx );
                if ( beta <= b_max and beta >= b_min ) gbr_thb[ nx ] += 1. * hbonded_t0[ bbnx ];
                if ( r <= rhbond_max and r >= rhbond_min  and beta <= betahbond_max and beta >= betahbond_min ) NHBt[ deltaTCF ] += 1. * hbonded_t0[ bbnx ];
                if ( r >  rhbond_max and r <= rhbond_maxT and beta <= betahbond_max and beta >= betahbond_min ) NTt[ deltaTCF ] += 1. * hbonded_t0[ bbnx ];
                if ( r <= rhbond_max and r >= rhbond_min  and beta >  betahbond_max and beta <= betahbond_maxR ) NRt[ deltaTCF ] += 1. * hbonded_t0[ bbnx ];

                // hydrogen 2
                bbnx = getbbnx( mol1, mol2, 2 );
                beta = bb[ bbnx ];
                bnx = get_bnx( beta );
                nx  = get_nx( rnx, bnx );
                if ( beta <= b_max and beta >= b_min ) gbr_thb[ nx ] += 1. * hbonded_t0[ bbnx ];
                if ( r <= rhbond_max and r >= rhbond_min  and beta <= betahbond_max and beta >= betahbond_min ) NHBt[ deltaTCF ] += 1. * hbonded_t0[ bbnx ];
                if ( r >  rhbond_max and r <= rhbond_maxT and beta <= betahbond_max and beta >= betahbond_min ) NTt[ deltaTCF ] += 1. * hbonded_t0[ bbnx ];
                if ( r <= rhbond_max and r >= rhbond_min  and beta >  betahbond_max and beta <= betahbond_maxR ) NRt[ deltaTCF ] += 1. * hbonded_t0[ bbnx ];
            }
        }
    }

    // normalization
    NHBt[ deltaTCF ]     /= 1.*nTCFsamples*nmol;
    NRt[ deltaTCF ]      /= 1.*nTCFsamples*nmol;
    NTt[ deltaTCF ]      /= 1.*nTCFsamples*nmol;

    // gbr_thb is normalized by gbr, but the number of samples must be the same for each
    for ( i = 0; i < npoints_r*npoints_b; i ++ ){
        // need to be careful here because if gbr[i] == 0, then gbr_thb will also be zero, but the division will be undefined
        if ( gbr[i] == 0 ) continue;
        else gbr_thb[i] /= 1.*gbr[i]*nTCFsamples/(1.*nsamples);
    }

    // may include option to write only certian ones that you want, but for now this is ok
    write_gbr_thb( deltaTCF );
}


void model::write_hbond_tcf()
{
    int deltaTCF;

    FILE *file = fopen(tcffname.c_str(),"w");
    fprintf( file, "#t (ns) hbond TCF\n");

    deltaTCF = 0;
    fprintf( file, "%g %g \n", deltaTCF*sampleEvery, hbondTCF[deltaTCF] );
    for ( int power = 0; power < 1000; power ++ ){
        if ( pow(10,power*.1) > deltaTCFMax ) break;
        deltaTCF = (int) round(pow(10,power*.1));
        if ( deltaTCF > deltaTCFMax ) break;
        fprintf( file, "%g %g \n", deltaTCF*sampleEvery, hbondTCF[deltaTCF] );
    }

    fclose( file );
}

void model::write_hbond_prt()
{
    int deltaTCF;

    FILE *file = fopen(prfname.c_str(),"w");
    fprintf( file, "#t (ns) hbond P_R\n");

    deltaTCF = 0;
    fprintf( file, "%g %g \n", deltaTCF*sampleEvery, NRt[deltaTCF] );
    for ( int power = 0; power < 1000; power ++ ){
        if ( pow(10,power*.1) > deltaTCFMax ) break;
        deltaTCF = (int) round(pow(10,power*.1));
        if ( deltaTCF > deltaTCFMax ) break;
        fprintf( file, "%g %g \n", deltaTCF*sampleEvery, NRt[deltaTCF]/NHBt[0] );
    }

    fclose( file );
}

void model::write_hbond_ptt()
{
    int deltaTCF;

    FILE *file = fopen(ptfname.c_str(),"w");
    fprintf( file, "#t (ns) hbond P_T\n");

    deltaTCF = 0;
    fprintf( file, "%g %g \n", deltaTCF*sampleEvery, NTt[deltaTCF] );
    for ( int power = 0; power < 1000; power ++ ){
        if ( pow(10,power*.1) > deltaTCFMax ) break;
        deltaTCF = (int) round(pow(10,power*.1));
        if ( deltaTCF > deltaTCFMax ) break;
        fprintf( file, "%g %g \n", deltaTCF*sampleEvery, NTt[deltaTCF]/NHBt[0] );
    }

    fclose( file );
}

void model::write_hbond_phbt()
{
    int deltaTCF;

    FILE *file = fopen(phbfname.c_str(),"w");
    fprintf( file, "#t (ns) hbond P_HB\n");

    deltaTCF = 0;
    fprintf( file, "%g %g \n", deltaTCF*sampleEvery, NHBt[deltaTCF] );
    for ( int power = 0; power < 1000; power ++ ){
        if ( pow(10,power*.1) > deltaTCFMax ) break;
        deltaTCF = (int) round(pow(10,power*.1));
        if ( deltaTCF > deltaTCFMax ) break;
        fprintf( file, "%g %g \n", deltaTCF*sampleEvery, NHBt[deltaTCF]/NHBt[0] );
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

    int currentSample, i, mol1, mol2, rnx, bnx, nx, rrnx, bbnx;
    float r, b1, b2, b;
    float ht0, NHB;
    float rho=0;

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

    // loop over trajectory
    for ( currentSample = 0; currentSample < reader.nsamples; currentSample ++ ){
        cout << "\rCurrent time: " << reader.gmxtime << setprecision(2) << fixed <<  " (ps)";
        cout.flush();

        // reset hbond array -- dont need to do this since each one is set explicitly
        memset( reader.hbonded, false, reader.nmol*reader.nmol*2*sizeof(bool));

        // calculate average density in molecules/nm
        rho += reader.nmol / (reader.box[0][0] * reader.box[1][1] * reader.box[2][2] ) / reader.nsamples ;
        
        // calculate H-bond distribution of angles and distances 
        // loop over all pairs of molecules
        for ( mol1 = 0; mol1 < reader.nmol; mol1 ++ ){
            for ( mol2 = 0; mol2 < reader.nmol; mol2 ++ ){

                if ( mol1 == mol2 ) continue;

                // Get OO distance
                rrnx = reader.getrrnx( mol1, mol2 );
                reader.rr[ rrnx ] = reader.get_r( mol1, mol2 );

                // Get Angle -- Hydrogen 1
                bbnx = reader.getbbnx( mol1, mol2, 1 );
                reader.bb[ bbnx ] = reader.get_b( mol1, mol2, 1 );
                if ( reader.bb[ bbnx ] < reader.b_max and reader.bb[ bbnx ] > reader.b_min ){
                    if ( reader.rr[ rrnx ] <= reader.r_max and reader.rr[ rrnx ] >= reader.r_min ){
                        rnx = reader.get_rnx( reader.rr[ rrnx ] );
                        bnx = reader.get_bnx( reader.bb[ bbnx ] );
                        nx  = reader.get_nx( rnx, bnx );
                        reader.gbr[ nx ] += 1.;
                    }
                }
                // determine if mol1 and mol2 O-OH pairs are hbonded for Hydrogen 1
                reader.hbonded[ bbnx ] = reader.is_hbond( reader.rr[ rrnx ], reader.bb[ bbnx ] );

                // Get Angle -- Hydrogen 2
                bbnx = reader.getbbnx( mol1, mol2, 2 );
                reader.bb[ bbnx ] = reader.get_b( mol1, mol2, 2 );
                if ( reader.bb[ bbnx ] < reader.b_max and reader.bb[ bbnx ] > reader.b_min ){
                    if ( reader.rr[ rrnx ] <= reader.r_max and reader.rr[ rrnx ] >= reader.r_min ){
                        rnx = reader.get_rnx( reader.rr[ rrnx ] );
                        bnx = reader.get_bnx( reader.bb[ bbnx ] );
                        nx  = reader.get_nx( rnx, bnx );
                        reader.gbr[ nx ] += 1.;
                    }
                }
                // determine if mol1 and mol2 O-OH pairs are hbonded for Hydrogen 2
                reader.hbonded[ bbnx ] = reader.is_hbond( reader.rr[ rrnx ], reader.bb[ bbnx ] );

            }
        }

        // write h-bond matrix, r matrix, and beta matrix to a scratch file
        reader.write_hbond(currentSample);

        // Advance to next frame if we need that frame
        if ( currentSample != reader.nsamples - 1 ) reader.search_for_sample( currentSample + 1 );
    } 

    // calculate the h-bond correlation function -- do this before 
    // normalization of gbr so you dont have to worry about that in the calculation of gbr(t|HB)
    int deltaSample = 0;
    cout << "\nNow calculating dynamics at t= " << deltaSample*reader.sampleEvery << " (ps)";
    cout.flush();
    reader.get_hbond_dynamics( deltaSample );
    // after calculating at time 0, just calculate on a log scale
    for ( int power = 0; power < 1000; power ++ ){
        if ( pow(10,power*.1) > reader.deltaTCFMax ) break;
            deltaSample = (int) round(pow(10,power*.1));

            if ( deltaSample > reader.deltaTCFMax ) break;
            cout << "\rNow calculating dynamics at t= " << deltaSample*reader.sampleEvery << " (ps)";
            cout.flush();
            reader.get_hbond_dynamics( deltaSample );
    }
    
    // remove hbond matrix scratch files
    reader.remove_hbond_files();
    
    // normalize hbond tcf
    for ( int i = reader.deltaTCFMax; i > 0; i -- ) reader.hbondTCF[i-1] /= reader.hbondTCF[0];

    NHB = 0;
    // normalize the probability distribution function, calculate the PMF, and the average number of Hbond acceptors
    for ( rnx = 0; rnx < reader.npoints_r; rnx ++ ){
        for ( bnx = 0; bnx < reader.npoints_b; bnx ++ ){
            nx = reader.get_nx( rnx, bnx );

            // normalize by the number of samples and the number of molecules
            reader.gbr[ nx ] /= reader.nsamples*reader.nmol*1.;

            // normalize by the "area" of each gridpoint -- take the center of each r, beta grid point.
            r = rnx * reader.dr + reader.r_min + reader.dr/2.;
            b = bnx * reader.db + reader.b_min + reader.db/2.;

            // integrate GBR over the h-bond region to find average number of hydrogen bonds
            if ( r <= reader.rhbond_max and r>= reader.rhbond_min and \
                    b <= reader.betahbond_max and b >= reader.betahbond_min ){
                NHB += reader.gbr[ nx ];
            }

            // normalize by the "area" of each gridpoint -- take the center of each r, beta grid point.
            // now gbr is independent of dr and dbeta
            reader.gbr[ nx ] /= 2. * PI * rho * r * r * sin( b * PI/180. ) * reader.dr * ( reader.db * PI / 180. );

            // note that the PMF is normalized by kT here
            reader.pmf[ nx ] = -1.*log(reader.gbr[nx]);
        }
    }

    cout << "\nThe density is " << rho << " molecules/nm3." << endl;
    cout << "The average number of hydrogen bonds is: " << NHB << setprecision(6) << endl;

    // write the probability distribution function and pmfs to a file
    reader.write_gbr();
    reader.write_pmf();
    reader.write_hbond_tcf();
    reader.write_hbond_prt();
    reader.write_hbond_ptt();
    reader.write_hbond_phbt();

    cout << endl << "DONE!" << endl;
}
