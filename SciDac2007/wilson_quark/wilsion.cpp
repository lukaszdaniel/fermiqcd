//#define SSE2
#define OSX
#include "fermiqcd.h"

#define X 0
#define Y 1
#define Z 2

void dump(mdp_field<float>& s,
	  int site_idx=0,
	  string filename="default.vtk",
	  bool bASCII=true){
  static const char *tempfile="tmp.vtk";
  char header[1024];
  char number[1024];
  char vectorHeader[512]; 
  int  LX=s.lattice().size(0), 
       LY=s.lattice().size(1), 
       LZ=s.lattice().size(2);
  FILE *file=NULL;
  
  snprintf(header, 1024,
	  "# vtk DataFile Version 2.0\n"
	  "Really cool data\n"
	  "%s\n"
	  "DATASET STRUCTURED_POINTS\n"
	  "DIMENSIONS %i %i %i\n"
	  "ORIGIN 0 0 0\n"
	  "SPACING 1 1 1\n"
	  "POINT_DATA %i\n"
	  "SCALARS scalar0 float 1\n"
  	  "LOOKUP_TABLE default\n",
	  ((bASCII)?"ASCII":"BINARY"),
	  LX+1,LY+1,LZ+1,
	  (LX+1)*(LY+1)*(LZ+1));

  remove(tempfile); //remove tmp.vtk file if exists
  file=fopen(tempfile, "wb");
  int sfd=fileno(file);
  cout << "saving... " << filename << " as fd=" << sfd << endl;
  fwrite(header, sizeof(char), strlen(header),file);
  site p(s.lattice());
  
  float fval=0.0F; 
  for(int i=0; i<LX+1; i++)
    for(int j=0; j<LY+1; j++)
      for(int k=0; k<LZ+1; k++){	
        p.set(i % LX,j % LY,k % LZ);

        fval=(float)s(p,site_idx);
	memset(number, 0, sizeof(number));
	snprintf(number, 1024,"%e\n",fval);
	fwrite(number, sizeof(char), strlen(number), file);

	/*
        fval=(float)s(p,1);
	memset(number, 0, sizeof(number));
	snprintf(number, 1024,"%e\n",fval);
	fwrite(number, sizeof(char), strlen(number), file);	
	*/
      }
 
  fclose(file);
  file=NULL;
  remove(filename.c_str());
  rename(tempfile, filename.c_str());
  return;
}


int main(int argc, char** argv) {
  define_base_matrices("FERMILAB");
  mdp.open_wormholes(argc,argv);
  int L[]={20,10,10,10};  
  mdp_lattice lattice(4,L);
  mdp_site x(lattice);
  gauge_field U(lattice,1);
  fermi_field psi(lattice,1);
  fermi_field phi(lattice,1);
  coefficients quark;
  quark["kappa"]=0.115;

  set_cold(U);
  psi=0;
  x.set(0,5,5,5);
  psi(x,0,0)=1; 

  mul_invQ(phi,psi,U,quark,1e-8);
  
  mdp_lattice space(3,L+1);
  mdp_field<float> s(space,9);
  mdp_site x3(space);
  char filename[128];
  for(int t=0; t<L[0]-1; t++) {
    forallsites(x3) {
      x.set(t % L[0],x3(0),x3(1),x3(2));
      s(x3,0)=real(phi(x,0,0));
      s(x3,1)=imag(phi(x,0,0));
      s(x3,2)=real(phi(x,1,0));
      s(x3,3)=imag(phi(x,1,0));
      s(x3,4)=real(phi(x,2,0));
      s(x3,5)=imag(phi(x,2,0));
      s(x3,6)=real(phi(x,3,0));
      s(x3,7)=imag(phi(x,3,0));
      s(x3,8)=sqrt(real(phi(x,0,0)*conj(phi(x,0,0))+phi(x,1,0)*conj(phi(x,1,0))+phi(x,2,0)*conj(phi(x,2,0))+phi(x,3,0)*conj(phi(x,3,0))));

    }
    snprintf(filename, 128,"wilson.s0.real.%.3i.vtk",t);
    dump(s,0,filename);
    snprintf(filename, 128,"wilson.s0.imag%.3i.vtk",t);
    dump(s,1,filename);
    snprintf(filename, 128,"wilson.s1.real.%.3i.vtk",t);
    dump(s,2,filename);
    snprintf(filename, 128,"wilson.s1.imag.%.3i.vtk",t);
    dump(s,3,filename);
    snprintf(filename, 128,"wilson.s2.real.%.3i.vtk",t);
    dump(s,4,filename);
    snprintf(filename, 128,"wilson.s2.imag.%.3i.vtk",t);
    dump(s,5,filename);
    snprintf(filename, 128,"wilson.s3.real.%.3i.vtk",t);
    dump(s,6,filename);
    snprintf(filename, 128,"wilson.s3.imag.%.3i.vtk",t);
    dump(s,7,filename);
    snprintf(filename, 128,"wilson.squared.%.3i.vtk",t);
    dump(s,8,filename);
  };

  mdp.close_wormholes();
  return 0;
}
