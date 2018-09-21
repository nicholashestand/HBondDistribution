#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <xdrfile.h>
#include <xdrfile_xtc.h>
#include <gmx_reader.h>
#include "HBondDistribution.h"

#define PI 3.14159265359 

using namespace std; 
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
model::model( string _inpf_ ) : gmx_reader::gmx_reader( _inpf_ )
// Default constructor
{

    // set userparams from input file
    for ( int i = 0; i < nuParams; i ++ )
    {
        if ( uParams[i] == "nsamples" )         nsamples        = stoi(uValues[i]);
        if ( uParams[i] == "sampleEvery" )      sampleEvery     = stof(uValues[i]);
        if ( uParams[i] == "beginTime" )        beginTime       = stof(uValues[i]);
        if ( uParams[i] == "deltaTCFMax" )      deltaTCFMax     = stoi(uValues[i]);
        if ( uParams[i] == "db" )               db              = stof(uValues[i]);
        if ( uParams[i] == "dr" )               dr              = stof(uValues[i]);
        if ( uParams[i] == "rmin" )             r_min           = stof(uValues[i]);
        if ( uParams[i] == "rmax" )             r_max           = stof(uValues[i]);
        if ( uParams[i] == "bmin" )             b_min           = stof(uValues[i]);
        if ( uParams[i] == "bmax" )             b_max           = stof(uValues[i]);
        if ( uParams[i] == "rhbond_max" )       rhbond_max      = stof(uValues[i]);
        if ( uParams[i] == "rhbond_min" )       rhbond_min      = stof(uValues[i]);
        if ( uParams[i] == "rhbondbreak_max" )  rhbondbreak_max = stof(uValues[i]);
        if ( uParams[i] == "betahbond_max")     betahbond_max   = stof(uValues[i]);
        if ( uParams[i] == "betahbond_min")     betahbond_min   = stof(uValues[i]);
        if ( uParams[i] == "betahbondbreak_max")betahbondbreak_max = stof(uValues[i]);
        if ( uParams[i] == "outfname" )         outfname        = uValues[i];
    }

    cout << "Set nsamples to: " << nsamples << endl;
    cout << "Set sampleEvery to: " << sampleEvery << endl;
    cout << "Set deltaTCFMax to: " << deltaTCFMax << endl;
    cout << "Set beginTime to: " << beginTime << endl;
    cout << "Set db to: "   << db << endl;
    cout << "Set dr to: "   << dr << endl;
    cout << "Set rmin to: " << r_min << endl;
    cout << "Set rmax to: " << r_max << endl;
    cout << "Set bmin to: " << b_min << endl;
    cout << "Set bmax to: " << b_max << endl;
    cout << "Set rhbond_max to: " << rhbond_max << endl;
    cout << "Set rhbond_min to: " << rhbond_min << endl;
    cout << "Set rhbondbreak_max" << rhbondbreak_max << endl;
    cout << "Set betahbond_max to: " << betahbond_max << endl;
    cout << "Set betahbond_min to: " << betahbond_min << endl;
    cout << "Set betahbondbreak_max to: " << betahbondbreak_max << endl;
    cout << "Set outfname to: " << outfname << endl;

    rhbondbreak_min = rhbond_max;
    betahbondbreak_min = betahbond_max;


    //  make sure deltaTCFMax is okay
    if ( deltaTCFMax > nframes*dt ){
        cout << "WARNING: deltaTCFMax > nframes*dt. Check input file" << endl;
        exit(1);
    }

    // determine number of gridpoints for r and beta
    npoints_r    = (int) round( (r_max-r_min)/dr + 1);           // total number of grid points
    npoints_b    = (int) round( (b_max-b_min)/db + 1);           // total number of grid points

    // determine the number of tcf points you will need so you can allocate enough memory
    // just care about points on a log scale and space them evenly
    ntcfpoints=0;
    for ( int tcfpoint = 0; tcfpoint < 1000; tcfpoint ++ )
    {
        float tcfdt = get_tcfdt( tcfpoint );
        if ( tcfdt > deltaTCFMax ) break;
        ntcfpoints ++;
    }
    
    // allocate arrays
    gbr         = new double[npoints_b*npoints_r]();
    pmf         = new double[npoints_b*npoints_r]();
    hbonded     = new bool[nmol*nmol*2]();
    hbondTCF    = new double[ntcfpoints]();
    NHBt        = new double[ntcfpoints]();
    NRt         = new double[ntcfpoints]();
    NTt         = new double[ntcfpoints]();
    gbr_thb     = new double[npoints_b*npoints_r*ntcfpoints]();
}

model::~model()
// Default Destructor
{
    delete [] gbr;
    delete [] pmf;
    delete [] hbonded;
    delete [] hbondTCF;
    delete [] gbr_thb;
}

float model::get_r( int mol1, int mol2 )
// Return the OO distance between two molecules
{
    float OOvec[3];
    int i;

    for ( i = 0; i < 3; i ++ ) OOvec[i] = x[ mol2 * natoms_mol ][i] \
                                        - x[ mol1 * natoms_mol ][i];
    minImage( OOvec );
    return mag3( OOvec );
}

float model::get_b( int mol1, int mol2, int whichH )
// Return the OH-O angle between two molecules
// OH is on molecule 1, O is on molecule 2
{
    float OOvec[3], OHvec[3];
    int i;
    float beta;
    float arg;

    if ( whichH < 1 or whichH > 2 ){
        cout << "Bad value for whichH: " << whichH << endl;
        exit(EXIT_FAILURE);
    }

    // OO vector
    for ( i = 0; i < 3; i ++ ) OOvec[i] = x[ mol2 * natoms_mol ][i] \
                                        - x[ mol1 * natoms_mol ][i];
    minImage( OOvec );

    // OH vector for the reference molecule
    for ( i = 0; i < 3; i ++ ) OHvec[i] = x[ mol1 * natoms_mol + whichH ][i] \
                                        - x[ mol1 * natoms_mol ][i];
    minImage( OHvec );

    // determine beta using trig
    arg   = dot3(OOvec,OHvec)/(mag3(OOvec)*mag3(OHvec));
    /* note acos is only defined when the argument is between zero and 1
     * -- this next line makes sure the program doesnt break due to numerical 
     * errors that cause the argument to go slightly over 1 or slightly below -1
     */
    arg   = max( arg, (float) -1. ); arg = min( arg, (float) 1. );
    if ( arg < -1 or arg > 1. ) cout << "WARNING: arg out of bounds" << endl;
    beta  = 180./PI*acos(arg);

    // return beta
    return  beta;
}

int model::get_rnx( float r )
// get grid index for a given value of r
{
    int rnx;

    rnx = round((r-r_min)/dr);
    // greater than or equal to because npoints_r includes 0
    if ( rnx < 0 or rnx >= npoints_r ) return -1;
    return rnx;
}

int model::get_bnx( float b )
// get grid indx for a given value of beta
{
    int bnx;
    
    bnx = round((b-b_min)/db);
    // greater than or equal to because npoints_b includes 0
    if ( bnx < 0 or bnx >= npoints_b ) return -1;
    return bnx;
}
                
int model::get_nx( int rnx, int bnx )
// get grid index for a rnx, bnx pair
{
    int nx;

    if ( rnx > npoints_r or rnx < 0) return -1;
    if ( bnx > npoints_b or bnx < 0) return -1;
    nx = bnx*npoints_r + rnx;
    return nx;
}

int model::get_nx_from_rb( float r, float b )
// get grid index for a r, beta pair of floats
{
    int rnx, bnx, nx;

    rnx = get_rnx( r );
    bnx = get_bnx( b );
    nx = get_nx( rnx, bnx ); 
    return nx;
}

float model::get_r_from_rnx( int rnx )
// get r from rnx
{
    float r;
    r = r_min + rnx*dr;
    return r;
}

float model::get_b_from_bnx( int bnx )
// get beta from bnx
{
    float b;
    b = b_min + bnx*db;
    return b;
}

int model::get_hbondnx( int mol1, int mol2, int h1 )
// and hbond index for the hbonded array
{
    if ( mol1 < 0 or mol1 > nmol ) return -1;
    if ( mol2 < 0 or mol2 > nmol ) return -1;
    if ( h1 < 1 or h1 > 2 ) return -1;
    return (mol1 * nmol + mol2)*2 + h1-1;
}

bool model::is_hbond( float r, float beta )
// determine if r and beta are hbonded or not
{
    // hbond condition is r<rhbond, beta<betahbond
    if ( r >= rhbond_max or r <= rhbond_min ) return false;
    if ( beta >= betahbond_max or beta <= betahbond_min ) return false;
    return true;
}

bool model::is_R_breakage( float r, float beta )
// determine if in rotational hbond breakage region
{
    if ( r >= rhbond_min and r <= rhbond_max and \
         beta > betahbondbreak_min and beta <= betahbondbreak_max ) return true;
    return false;
}

bool model::is_T_breakage( float r, float beta )
// determine if in rotational hbond breakage region
{
    if ( r >  rhbondbreak_min and r <= rhbondbreak_max and \
         beta >= betahbond_min and beta <= betahbond_max ) return true;
    return false;
}

void model::determine_hbonds()
// determine hbonds for the current frame
{
    int mol1, mol2, nx, i;
    float r, b;

    // initialize the hbonded array to false
    for ( nx = 0; nx < nmol*nmol*2 ; nx ++ ) hbonded[ nx ] = false;

    for ( mol1 = 0; mol1 < nmol; mol1 ++ ){
        for ( mol2 = 0; mol2 < nmol; mol2 ++ ){

            if ( mol1 == mol2 ) continue;

            r = get_r( mol1, mol2 );
            // speed it up a little bit -- if we dont care just skip it
            if ( r > r_max ) continue;
            
            // check both H for being hbonded
            for ( i = 1; i < 3; i ++ ){
                b = get_b( mol1, mol2, i);
                nx = get_hbondnx( mol1, mol2, i );
                hbonded[nx] = is_hbond( r, b );
            }
        }
    }
}

void model::determine_gbr()
// determine gbr for the current frame
{
    int mol1, mol2, nx, i, rnx, bnx;
    float r, b;

    for ( mol1 = 0; mol1 < nmol; mol1 ++ ){
        for ( mol2 = 0; mol2 < nmol; mol2 ++ ){
            if ( mol1 == mol2 ) continue;
            
            r = get_r( mol1, mol2 );
            // speed it up a little bit -- if we dont care just skip it
            if ( r > r_max ) continue;
            
            for ( i = 1; i < 3; i ++ ){
                b = get_b( mol1, mol2, i);
                nx = get_nx_from_rb( r, b );
                if ( nx != -1 ) gbr[ nx ] += 1.;
                // -1 returned by get_nx_from_rb if r and b are out of bounds
            }
        }
    }
}

float model::get_tcfdt( int tcfpoint )
// return the time offset for a give tcf point index
{
    float tcfdt;

    // let the first ten points be at the resolution of dt to avoid overlapping
    if ( tcfpoint < 10 ) tcfdt = tcfpoint*dt;
    // afterwards, increase tcfdt on a long scale with 10 points per power of 10
    else tcfdt = round(pow(10,(tcfpoint)*0.1))*dt;
    return tcfdt;
}


void model::calculate_dynamicTCFs( int tcfpoint )
// calculate all of the dynamical quantities of interest
{

    int mol1, mol2, nx, hbnx;
    float r, b, h;

    // loop over all pairs of molecules
    for ( mol1 = 0; mol1 < nmol; mol1 ++ ){
        for ( mol2 = 0; mol2 < nmol; mol2 ++ ){

            if ( mol1 == mol2 ) continue;

            r = get_r( mol1, mol2 );
            // speed it up a little bit -- if we dont care just skip it
            if ( r > r_max ) continue;

            // hbond correlation function
            for ( h = 1; h < 3; h ++ ){
                hbnx = get_hbondnx( mol1, mol2, h );
                // speed it up a little more -- below will only contribute if hbonded[ hbnx ] is true
                if ( hbonded[ hbnx ] ){
                    b = get_b( mol1, mol2, h);

                    // correlation functions -- note these are accumulated for each tcfpoint
                    hbondTCF[ tcfpoint ] += hbonded[ hbnx ] * is_hbond( r, b );
                    NHBt[ tcfpoint ]     += hbonded[ hbnx ] * is_hbond( r, b );
                    NRt[ tcfpoint ]      += hbonded[ hbnx ] * is_R_breakage( r, b );
                    NTt[ tcfpoint ]      += hbonded[ hbnx ] * is_T_breakage( r, b );

                    // gbr(thb)
                    nx = get_nx_from_rb( r, b );
                    // note that the index here is a little different because at each time
                    // step there are npoints_r*nponts_b possible entries
                    if ( nx != -1 ) gbr_thb[ tcfpoint*npoints_r*npoints_b + nx ] += hbonded[ hbnx ];
                }
            }

        }
    }
}

void model::write_gbr()
// write average gbr to file
{
    int rnx, bnx, nx;
    float r, b;

    string fname = outfname+"-gbr.dat";
    FILE *file = fopen(fname.c_str(),"w");
    fprintf( file, "#R (nm) Beta (deg) g(R,Beta)\n");

    for ( rnx = 0; rnx < npoints_r; rnx ++ ){
        r = get_r_from_rnx( rnx );
        for ( bnx = 0; bnx < npoints_b; bnx ++ ){
            nx = get_nx( rnx, bnx );
            b = get_b_from_bnx( bnx );
            fprintf( file, "%g %g %g\n", r, b, gbr[ nx ]);
        }
    }
    fclose( file );
}

void model::write_pmf()
// write average pmf to file
{
    int rnx, bnx, nx;
    float r, b;

    string fname = outfname+"-pmf.dat";
    FILE *file = fopen(fname.c_str(),"w");
    fprintf( file, "#R (nm) Beta (deg) PMF(R,Beta)\n");

    for ( rnx = 0; rnx < npoints_r; rnx ++ ){
        for ( bnx = 0; bnx < npoints_b; bnx ++ ){
            nx = get_nx( rnx, bnx );
            r = get_r_from_rnx( rnx );
            b = get_b_from_bnx( bnx );
            fprintf( file, "%g %g %g\n", r, b, pmf[ nx ]);
        }
    }
    fclose( file );
}


void model::write_dynamicTCFs()
{
    int tcfpoint, rnx, bnx, nx;
    FILE *file;
    string fname;
    float r, b, tcfdt;

    // hbond tcf
    fname = outfname+"-hbtcf.dat";
    file = fopen(fname.c_str(),"w");
    fprintf( file, "#t (ps) hbond TCF\n");
    for ( tcfpoint = 0; tcfpoint < ntcfpoints; tcfpoint ++ ){
        tcfdt = get_tcfdt( tcfpoint );
        fprintf( file, "%g %g \n", tcfdt, hbondTCF[ tcfpoint ]/hbondTCF[0] );
    }
    fclose( file );

    // PRt
    fname = outfname+"-nrt.dat";
    file = fopen(fname.c_str(),"w");
    fprintf( file, "#t (ps) P_R\n");
    for ( tcfpoint = 0; tcfpoint < ntcfpoints; tcfpoint ++ ){
        tcfdt = get_tcfdt( tcfpoint );
        fprintf( file, "%g %g \n", tcfdt, NRt[ tcfpoint ]/NHBt[0] );
    }
    fclose( file );

    // PTt
    fname = outfname+"-ntt.dat";
    file = fopen(fname.c_str(),"w");
    fprintf( file, "#t (ps) P_T\n");
    for ( tcfpoint = 0; tcfpoint < ntcfpoints; tcfpoint ++ ){
        tcfdt = get_tcfdt( tcfpoint );
        fprintf( file, "%g %g \n", tcfdt, NTt[ tcfpoint ]/NHBt[0] );
    }
    fclose( file );

    // gbrthb
    for ( tcfpoint = 0; tcfpoint < ntcfpoints; tcfpoint ++ ){
        fname = outfname+"-gbrthb-"+to_string(tcfpoint)+".dat";
        file = fopen(fname.c_str(),"w");
        tcfdt = get_tcfdt( tcfpoint );
        fprintf( file, "#Time: %g (ps) \n", tcfdt );
        fprintf( file, "#R (nm) Beta (deg) g(R,Beta,t)\n");
        for ( rnx = 0; rnx < npoints_r; rnx ++ ){
            for ( bnx = 0; bnx < npoints_b; bnx ++ ){
                nx = get_nx( rnx, bnx );
                b = get_b_from_bnx( bnx );
                r = get_r_from_rnx( rnx );
                fprintf( file, "%g %g %g\n", r, b, gbr_thb[ tcfpoint*npoints_r*npoints_b + nx ]);
            }
        }
        fclose(file);
    }
}


// ************************************************************************ 
int main( int argc, char* argv[] )
{

    int currentSample, i, mol1, mol2, rnx, bnx, nx, rrnx, bbnx;
    int frameno, tcfpoint;
    float r, b;
    float rho=0, NHB;
    float currentTime, tcfTime, tcfdt;

    // Check program input
    if ( argc != 2 ){
        printf("Program expects only one argument, which is the name of \n\
                an input file containing the details of the analysis.\nAborting...\n");
        exit(EXIT_FAILURE);
    }

    // get filename for parameters
    string inpf(argv[1]);

    // initialize class 
    model reader( inpf );

    // take samples for hbond dynamics
    cout << endl << "************************************************************************" << endl;
    for ( currentSample = 0; currentSample < reader.nsamples; currentSample ++ ){

        // read frame at t0
        currentTime = currentSample * reader.sampleEvery + reader.beginTime;
        cout << "\rCurrent time: " << currentTime << setprecision(2) << fixed <<  " (ps)";
        cout.flush();
        frameno = reader.get_frame_number( currentTime );
        if ( frameno == -1 ){
            cout << "Warning:: get_frame_number returned " << frameno << " when searching for frame at " \
                << currentTime << " (ps).\nCheck input file that nsamples is not too big!" << endl;
            exit(EXIT_FAILURE);
        }
        reader.read_frame( frameno );
        if ( reader.checktime( currentTime ) == false ){
            cout << "Warning:: checktime failed. gmxtime is: " << reader.gmxtime << endl;
            cout << "Aborting." << endl;
            exit(EXIT_FAILURE);
        }

        // calculate average density in molecules/nm
        rho += reader.nmol / (reader.box[0][0] * reader.box[1][1] * reader.box[2][2] ) / reader.nsamples ;

        // determine hbonds and gbr at t0
        reader.determine_hbonds(); 
        reader.determine_gbr();

        // loop over offsets to find time correlations
        for ( tcfpoint = 0; tcfpoint < reader.ntcfpoints; tcfpoint ++ )
        {
            tcfdt = reader.get_tcfdt( tcfpoint ); 
            tcfTime = currentTime + tcfdt;
            frameno = reader.get_frame_number( tcfTime );
            reader.read_frame( frameno );
            reader.calculate_dynamicTCFs( tcfpoint );
        }
    }   

    // normalize the gbr_thb functions
    for ( tcfpoint = 0; tcfpoint < reader.ntcfpoints; tcfpoint ++ ){
        for ( i = 0; i < reader.npoints_r*reader.npoints_b; i ++ ){
            reader.gbr_thb[ tcfpoint*reader.npoints_r*reader.npoints_b + i ] /= reader.gbr[i];
        }
    }

    NHB = 0;
    // normalize gbr and calculate the pmf
    for ( rnx = 0; rnx < reader.npoints_r; rnx ++ ){
        for ( bnx = 0; bnx < reader.npoints_b; bnx ++ ){
            nx = reader.get_nx( rnx, bnx );
            r  = reader.get_r_from_rnx( rnx );
            b  = reader.get_b_from_bnx( bnx );

            // normalize by the number of samples and the number of molecules
            reader.gbr[ nx ] /= reader.nsamples*reader.nmol*1.;

            // integrate GBR over the h-bond region to find average number of hydrogen bonds
            if ( r <= reader.rhbond_max and r>= reader.rhbond_min and \
                    b <= reader.betahbond_max and b >= reader.betahbond_min ){
                NHB += reader.gbr[ nx ];
            }

            /* normalize by the "area" of each gridpoint -- 
             * take the center of each r, beta grid point.
             * now gbr is independent of dr and dbeta */
            reader.gbr[ nx ] /= 2. * PI * rho * r * r * sin( b * PI/180. ) * \
                                reader.dr * ( reader.db * PI / 180. );

            // note that the PMF is normalized by kT here
            reader.pmf[ nx ] = -1.*log(reader.gbr[nx]);
        }
    }

    cout << "\nThe density is " << rho << " molecules/nm3." << endl;
    cout << "The average number of hydrogen bonds is: " << NHB << setprecision(6) << endl;

    // write the probability distribution function and pmfs to a file
    reader.write_gbr();
    reader.write_pmf();
    reader.write_dynamicTCFs();

    cout << endl << "DONE!" << endl;
}
