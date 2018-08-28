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
        if ( uParams[i] == "gbrname" )  gbrfname = uValues[i];
        if ( uParams[i] == "pmfname" )  pmffname = uValues[i];

    }

    cout << "Set db to: "   << db << endl;
    cout << "Set dr to: "   << dr << endl;
    cout << "Set rmin to: " << r_min << endl;
    cout << "Set rmax to: " << r_max << endl;
    cout << "Set bmin to: " << b_min << endl;
    cout << "Set bmax to: " << b_max << endl;
    cout << "Set gbrname to: " << gbrfname << endl;
    cout << "Set pmfname to: " << pmffname << endl;

    // allocate space for nearest neighbor classifications
    npoints_r    = (int) round( (r_max-r_min)/dr + 1);           // total number of grid points
    npoints_b    = (int) round( (b_max-b_min)/db + 1);           // total number of grid points

    // the array for the H-bond distribution function
    gbr    = new double[npoints_b*npoints_r]();
    pmf    = new double[npoints_b*npoints_r]();
}

// Default Destructor
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
model::~model()
{
    delete [] gbr;
    delete [] pmf;
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

    // return the minimum beta angle
    arg   = dot3(OOvec,OHvec)/(mag3(OOvec)*mag3(OHvec));
    arg   = max( arg, (float) -1. ); arg = min( arg, (float) 1. );
    if ( arg < -1 or arg > 1. ) cout << "WARNING: arg out of bounds" << endl;

    beta  = 180./PI*acos(arg);
    // i have to do this because acos has trouble with when the argument is 1, it returns nan
    //if ( isnan(beta1) ) beta1 = 0.;

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


int main( int argc, char* argv[] )
{

    int currentSample, i, mol1, mol2, rnx, bnx, nx;
    long int kount;
    float r, b1, b2;

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

        // calculate H-bond distribution of angles and distances 
        // loop over all pairs of molecules
        for ( mol1 = 0; mol1 < reader.nmol; mol1 ++ ){
            for ( mol2 = 0; mol2 < reader.nmol; mol2 ++ ){
                if ( mol1 == mol2 ) continue;

                // increment counter -- must divide at the end to get a probability -- each water has 2 OH angles
                kount += 2;

                // Get OO distance
                r = reader.get_r( mol1, mol2 );
                if ( r > reader.r_max or r < reader.r_min ) continue;

                // Get Angle -- Hydrogen 1
                b1 = reader.get_b( mol1, mol2, 1 );
                if ( b1 < reader.b_max and b1 > reader.b_min ){
                    rnx = reader.get_rnx( r );
                    bnx = reader.get_bnx( b1 );
                    nx  = reader.get_nx( rnx, bnx );
                    reader.gbr[ nx ] += 1.;
                }

                // Get Angle -- Hydrogen 2
                b2 = reader.get_b( mol1, mol2, 2 );
                if ( b2 < reader.b_max and b2 > reader.b_min ){
                    rnx = reader.get_rnx( r );
                    bnx = reader.get_bnx( b2 );
                    nx  = reader.get_nx( rnx, bnx );
                    reader.gbr[ nx ] += 1.;
                }
            }
        }
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

    // write the probability distribution function and pmfs to a file
    reader.write_gbr();
    reader.write_pmf();

    cout << endl << "DONE!" << endl;
}
