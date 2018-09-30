#include<stdio.h>
#include<math.h>

#define X 64
#define Y X
#define MR 0.125

#define LX 3
#define LY LX
#define T 0.5
#define dX (((double)LX)/X)
#define dY (((double)LY)/Y)

#define EX 1

double f(double u);

double dfdu(double u);

double g(double u);

double dgdu(double u);

void KT(double u[X+3][Y+3], double Up[X+3][Y+3], double Um[X+3][Y+3], double Vp[X+3][Y+3], double Vm[X+3][Y+3], double phi);

void ExTest();

void Ex();

double minmod(double a, double b, double c);

double sgn(double a);

double minabs(double a, double b);

double max(double a, double b);

double min(double a, double b);

void printResults(double u[X+3][Y+3], char filename[]);

int main(){
  switch(EX){
    case 1: ExTest(); break;
    case 2: Ex(); break;
  }
  return 0;
}

void KT(double u[X+3][Y+3], double Up[X+3][Y+3], double Um[X+3][Y+3], double Vp[X+3][Y+3], double Vm[X+3][Y+3], double phi){
  double up[X+3][X+3];
  double ux[X+3][X+3], uy[X+3][X+3];
  double uxpp[X+3][X+3], uxpm[X+3][X+3], uxmp[X+3][X+3], uxmm[X+3][X+3];
  double uypp[X+3][X+3], uypm[X+3][X+3], uymp[X+3][X+3], uymm[X+3][X+3];
  double axp[X+3][X+3], axm[X+3][X+3];
  double ayp[X+3][X+3], aym[X+3][X+3];
  double Hxplus[X+3][X+3], Hxminus[X+3][X+3];
  double Hyplus[X+3][X+3], Hyminus[X+3][X+3];
  double du[X+3][X+3];

  double dT = MR*(dX);
  double lbd = MR;
  double theta = 2;
  int t, j, k;

  char filename[30];

  //time :: Euler
  //u_j(t+dT) = u_j(t) + ((dU/dt)_j)*dT
  double elp = 0;
  for(t = 0; elp < T; t++){ // [pretty done]
    //ux,uy vector calculation [done]
    for(k = 1; k <= X+1; k++){
	    for(j = 1; j <= X+1; j++){
	      ux[j][k] = minmod(theta*(u[j][k]-u[j-1][k])/(dX), (u[j+1][k]-u[j-1][k])/(2*dX), theta*(u[j+1][k]-u[j][k])/(dX));
	    }
	  }
    for(j = 1; j <= X+1; j++){
	    for(k = 1; k <= X+1; k++){
	      uy[j][k] = minmod(theta*(u[j][k]-u[j][k-1])/(dX), (u[j][k+1]-u[j][k-1])/(2*dX), theta*(u[j][k+1]-u[j][k])/(dX));
	    }
	  }

    //uxxx and uyyy calculation [done]
    for(k = 1; k <= X+1; k++){
      for(j = 1; j <= X+1; j++){
        uxpp[j][k] = u[j+1][k]-((dX)/2)*(ux[j+1][k]);
        uxmp[j][k] = u[j][k]+((dX)/2)*(ux[j][k]);
        uxpm[j][k] = u[j][k]-((dX)/2)*(ux[j][k]);
        uxmm[j][k] = u[j-1][k]+((dX)/2)*(ux[j-1][k]);
      }
    }
    for(j = 1; j <= X+1; j++){
      for(k = 1; k <= X+1; k++){
        uypp[j][k] = u[j][k+1]-((dY)/2)*(uy[j][k+1]);
        uymp[j][k] = u[j][k]+((dY)/2)*(uy[j][k]);
        uypm[j][k] = u[j][k]-((dY)/2)*(uy[j][k]);
        uymm[j][k] = u[j][k-1]+((dY)/2)*(uy[j][k-1]);
      }
    }

    //a calculation [done]
    for(k = 1; k <= X+1; k++){
      for(j = 1; j <= X+1; j++){
        axp[j][k] = max(fabs((Up[j][k]/phi)*dfdu(uxpp[j][k])), fabs((Up[j][k]/phi)*dfdu(uxmp[j][k])));
        axm[j][k] = max(fabs((Um[j][k]/phi)*dfdu(uxpm[j][k])), fabs((Um[j][k]/phi)*dfdu(uxmm[j][k])));
        ayp[j][k] = max(fabs((Vp[j][k]/phi)*dgdu(uypp[j][k])), fabs((Vp[j][k]/phi)*dgdu(uymp[j][k])));
        aym[j][k] = max(fabs((Vm[j][k]/phi)*dgdu(uypm[j][k])), fabs((Vm[j][k]/phi)*dgdu(uymm[j][k])));
      }
    }

    //H calculation (xminus,xplus,yminus,yplus) [done]
    for(k = 1; k <= X+1; k++){
      for(j = 1; j <= X+1; j++){
        Hxplus[j][k] = (Up[j][k]*f(uxpp[j][k])+Up[j][k]*f(uxmp[j][k]))/(2*phi) - (axp[j][k]/2)*(uxpp[j][k]-uxmp[j][k]);
        Hxminus[j][k] = (Um[j][k]*f(uxpm[j][k])+Um[j][k]*f(uxmm[j][k]))/(2*phi) - (axm[j][k]/2)*(uxpm[j][k]-uxmm[j][k]);
        Hyplus[j][k] = (Vp[j][k]*g(uypp[j][k])+Vp[j][k]*g(uymp[j][k]))/(2*phi) - (ayp[j][k]/2)*(uypp[j][k]-uymp[j][k]);
        Hyminus[j][k] = (Vm[j][k]*g(uypm[j][k])+Vm[j][k]*g(uymm[j][k]))/(2*phi) - (aym[j][k]/2)*(uypm[j][k]-uymm[j][k]);
      }
    }

    //du calculation [done]
    for(k = 1; k <= X+1; k++){
      for(j = 1; j <= X+1; j++){
        du[j][k] = (Hxminus[j][k]-Hxplus[j][k])/(dX) + (Hyminus[j][k]-Hyplus[j][k])/(dY);
      }
    }

    //u calculation (Euler) [done]
    for(k = 1; k <= X+1; k++){
      for(j = 1; j <= X+1; j++){
        u[j][k] = u[j][k] + du[j][k]*(dT);
      }
      //u[0][k] = u[1][k];
      //u[X+2][k] = u[X+1][k];
      /*instead if doing this, fill boundary cells with
      some kind of default-value (1 for max saturation)*/
    }

    sprintf(filename, "temp/%d.txt", t);
    printResults(u, filename);

    elp += dT; //increase elapsed time
  }

}

void ExTest(){

  double u[X+3][X+3];
  double Up[X+3][X+3], Um[X+3][X+3], Vp[X+3][X+3], Vm[X+3][X+3];
  double phi = 1;
  int j, k;

  //initial conditions
  for(k = 1; k <= X+1; k++){
    for(j = 1; j <= X+1; j++){
      if(pow(j*dX-1.0, 2)+pow(k*dY-1.0, 2) <= pow(0.4, 2)){
        u[j][k] = 1;
      } else if(pow(j*dX-2.0, 2)+pow(k*dY-2.0, 2) <= pow(0.4, 2)){
        u[j][k] = -1;
      }
    }
  }
  //-------------------

  //boundary conditions
  for(k = 0; k <= Y+2; k++){
    u[0][k] = 0;
    u[X+2][k] = 0;
  }
  for(j = 0; j <= X+2; j++){
    u[j][0] = 0;
    u[j][X+2] = 0;
  }
  //-------------------

  //U,V values
  for(k = 0; k <= X+2; k++){
    for(j = 0; j <= X+2; j++){
      Up[j][k] = 1;
      Um[j][k] = 1;
      Vp[j][k] = 1;
      Vm[j][k] = 1;
    }
  }
  //-------------------

  KT(u, Up, Um, Vp, Vm, phi);

  printResults(u, "ExTestKT.txt");
}


void Ex(){

  double u[X+3];
  int j;

  //load vel. values


  //-------------------
  //set U,V value

  //-------------------
  //set boundary conditions
  for(j = 1; j <= X+1; j++){
    u[j] = sin((j-1)*dX);
  }
  u[0] = u[1]; //ghost
  u[X+2] = u[X+1]; //ghost
  //-------------------

  //KT(u);
}

//########### f(u) ###########//
double f(double u){
  switch(EX){
    case 1: return pow(u, 2); break;
    case 2: return 0; break;
  }
  return 0;
}
//###########################//
//########## df/du ##########//
double dfdu(double u){
  switch(EX){
    case 1: return 2*u; break;
    case 2: return 0; break;
  }
  return 0;
}
//########### g(u) ###########//
double g(double u){
  switch(EX){
    case 1: return pow(u, 2); break;
    case 2: return 0; break;
  }
  return 0;
}
//###########################//
//########## dg/du ##########//
double dgdu(double u){
  switch(EX){
    case 1: return 2*u; break;
    case 2: return 0; break;
  }
  return 0;
}
//###########################//

//###########################//
//           UTILS           //
//###########################//

double minmod(double a, double b, double c){
  if(a < 0 && b < 0 && c < 0){
    return max(max(a, b), c);
  } else if(a > 0 && b > 0 && c > 0){
    return min(min(a, b), c);
  } else {
    return 0;
  }
}

double sgn(double a){
  return (0 < a) - (a < 0);
}

double minabs(double a, double b){
  if(fabs(a) < fabs(b)) return fabs(a);
  else return fabs(b);
}

double max(double a, double b){
  if(a > b) return a;
  else return b;
}

double min(double a, double b){
  if(a < b) return a;
  else return b;
}

void printResults(double u[X+3][Y+3], char filename[]){
  int j, k;
  FILE* f = fopen(filename, "w");
  for(k = 1; k <= Y+1; k++){
    for(j = 1; j <= X+1; j++){
      fprintf(f, "%lf %lf %lf\n", (j-1)*dX, (k-1)*dY, u[j][k]);
    }
  }
  fclose(f);
}

//###########################//
//           END ;)          //
//###########################//
