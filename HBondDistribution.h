// create a new subclase that inherits the gmx reader class
#include <string.h>
#include <gmx_reader.h>

#ifndef HBondDistribution_H
#define HBondDistribution_H
class model: public gmx_reader
{
    public:
        // class variables
        float db=1;                                             // grid spacing in degrees
        float dr=0.1;                                           // grid spacing in nm
        float r_min=0., r_max=0.6;                              // max and min grid points in nm
        float b_min=0., b_max=90;                               // max and min grid points in degree
        int npoints_r;                                          // total number of grid points
        int npoints_b;                                          // total number of grid points
        float rhbond;
        float betahbond;
        double *gbr, *pmf;                                      // the probability distribution function
        bool *hbonded,*hbonded_t0,*hbonded_t;                   // boolean array to keep track of hbonds between pairs of molecules
        float *hbondTCF;                                        // hydrogenbond correlation function
        string gbrfname="gbr.dat", pmffname="pmf.dat";          // filenames for outputs
        string tcffname="hbondtcf.dat";                         // more filenames
        int deltaTCFMax = 10;
        float *rr, *bb1, *bb2;                                  // r, beta1, and beta 2 for each molecule
        double *gbr_thb;                                        // instantaneous distribution at time t

    
        // Default constructor and destructor
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        model( string _inpf_ );
        ~model();

        // Nearest neighbor functions
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        float get_r( int mol1, int mol2 );
        float get_b( int mol1, int mol2, int whichH );
        int   get_rnx( float r );
        int   get_bnx( float b );
        int   get_nx( int rnx, int bnx );
        void  write_gbr( );
        void  write_pmf( );
        void  write_hbond_tcf( );
        bool  is_hbond( float r, float beta1, float beta2 );
        void  write_hbond( int currentSample );
        void  remove_hbond_files();
        void  read_hbond_t0( int currentSample );
        void  read_hbond_t( int currentSample );
        float get_hbond_TCF( int deltaSample );
        int   getarraynx( int mol1, int mol2 );

};
#endif
