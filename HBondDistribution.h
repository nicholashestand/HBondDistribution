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
        double *gbr, *pmf;                                      // the probability distribution funciton
        string gbrfname="gbr.dat", pmffname="pmf.dat";          // filenames for outputs

    
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

};
#endif
