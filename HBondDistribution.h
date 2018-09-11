// create a new subclase that inherits the gmx reader class
#include <string.h>
#include <gmx_reader.h>

#ifndef HBondDistribution_H
#define HBondDistribution_H
class model: public gmx_reader
{
    public:
        // class variables
        string outfname="hbond";
        float db=1, dr=0.1, r_min=0, r_max=0.6, b_min=0, b_max=90;  // grid spacing, and bounds
        int npoints_r, npoints_b, ntcfpoints;                       // total number of grid points
        float rhbond_min, rhbond_max;                               // parameters to define hbonds
        float betahbond_min, betahbond_max, betahbondbreak_max;
        float betahbondbreak_min, rhbondbreak_min, rhbondbreak_max;
        float deltaTCFMax;
        double *gbr, *pmf;                                          // the probability distribution function
        bool   *hbonded_t0,*hbonded;                                // boolean array to keep track of hbonds between OHO pairs
        double *hbondTCF, *NHBt, *NRt, *NTt;                        // correlation functions
        int nsamples, sampleEvery;                                  // sampling variables
        int nTCFpoints;                                             // number of TCF points -- calculated internally
        double *gbr_thb;


    
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
        int   get_nx_from_rb( float r, float b );
        float get_r_from_rnx( int rnx );
        float get_b_from_bnx( int bnx );
        int   get_hbondnx( int mol1, int mol2, int h );
        bool  is_hbond( float r, float beta);
        bool  is_R_breakage( float r, float beta);
        bool  is_T_breakage( float r, float beta);
        void  determine_hbonds();
        void  determine_gbr();
        float get_tcfdt( int tcfpoint );
        void  calculate_dynamicTCFs( int tcfpoint );
        void  write_gbr( );
        void  write_pmf( );
        void  write_dynamicTCFs();
};
#endif
