#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main()

{

	FILE *ifp;
	FILE *ofp;

	ifp = fopen("dump.lammpstrj.physicalwvir","r");
	ofp = fopen("cg.lammpstrj.physicalwvir","r");

	double x[10];
	double y[10];
	double z[10];

	double fx[10];
	double fy[10];
	double fz[10];

	double xcg;
	double ycg;
	double zcg;

	double fxcg;
	double fycg;
	double fzcg;

	double disx;
	double disy;
	double disz;
	
	double nbox;
	double box;
	
	int t;
	int i;
	int j;
	int k;
	int id;
	int type;
	int ts;

	char buffer[100];

	for(t=0;t<10;t++){

		printf("%d\n",t);

		fgets(buffer,100,ifp);
		fprintf(ofp,"ITEM: TIMESTEP\n",buffer);
		
		fscanf(ifp,"%d\n",&ts);
		fprintf(ofp,"%d\n",ts);

		fgets(buffer,100,ifp);
		fprintf(ofp,"ITEM: NUMBER OF ATOMS\n");

		fgets(buffer,100,ifp);
		fprintf(ofp,"250\n");

		fgets(buffer,100,ifp);
		fprintf(ofp,"ITEM: BOX BOUNDS pp pp pp\n");

		fscanf(ifp,"%lf %lf\n",&nbox,&box);
		fprintf(ofp,"0 31.66\n");

		fgets(buffer,100,ifp);
		fprintf(ofp,"0 31.66\n");

		fgets(buffer,100,ifp);
		fprintf(ofp,"0 31.66\n");

		fgets(buffer,100,ifp);
		fprintf(ofp,"ITEM: ATOMS id type x y z fx fy fz \n");
		printf("%f\n",box);

		for(i=0;i<250;i++){

			//printf("%d\n",i);

			for(j=0;j<=9;j++){
				fscanf(ifp,"%d %d %lf %lf %lf %lf %lf %lf\n",&id,&type,&x[j],&y[j],&z[j],&fx[j],&fy[j],&fz[j]);
			}
			
			//printf("check two\n");
			for(j=1;j<=9;j++){

				disx = x[j]-x[0];
				disy = y[j]-y[0];
				disz = z[j]-z[0];

				if(disx > box/2.0){
				x[j] = x[j] - box;
				}
				if(disx < -box/2.0){
				x[j] = x[j] + box;
				}
				if(disy > box/2.0){
				y[j] = y[j] - box;
				}
				if(disy < -box/2.0){
				y[j] = y[j] + box;
				}	
				if(disz > box/2.0){
				z[j] = z[j] - box;
				}
				if(disz < -box/2.0){
				z[j] = z[j] + box;
				}
			}			

			//printf("check three\n");

			//calculate COM
			xcg = (15.999400*x[0]+15.999400*x[1]+15.999400*x[2]+12.010700*x[3]+12.010700*x[4]+12.010700*x[5]+1.007940*x[6]+1.007940*x[7]+1.007940*x[8]+1.007940*x[9])/88.06206;
			ycg = (15.999400*y[0]+15.999400*y[1]+15.999400*y[2]+12.010700*y[3]+12.010700*y[4]+12.010700*y[5]+1.007940*y[6]+1.007940*y[7]+1.007940*y[8]+1.007940*y[9])/88.06206;
			zcg = (15.999400*z[0]+15.999400*z[1]+15.999400*z[2]+12.010700*z[3]+12.010700*z[4]+12.010700*z[5]+1.007940*z[6]+1.007940*z[7]+1.007940*z[8]+1.007940*z[9])/88.06206;

			//printf("check four\n");


			if(xcg<0.0){
			xcg = xcg + box;
			}
			if(xcg>box){
			xcg = xcg - box;
			}
			if(ycg<0.0){
			ycg = ycg + box;
			}
			if(ycg>box){
			ycg = ycg - box;
			}
			if(zcg<0.0){
			zcg = zcg + box;
			}
			if(zcg>box){
			zcg= zcg -box;
			}
		

			fxcg = 0.0;
			fycg = 0.0;
			fzcg = 0.0;
			
			for(j=0;j<=9;j++){
				fxcg= fxcg+fx[j];
				fycg= fycg+fy[j];
				fzcg= fzcg+fz[j];		
			}

			//printf("check five\n");

			fprintf(ofp,"%d %d %.6g %.6g %.6g %.6g %.6g %.6g\n",i+1,1,xcg,ycg,zcg,fxcg,fycg,fzcg);

			//printf("check six\n");
		}


	}
	return 0;
}
