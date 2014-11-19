#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main()

{

	FILE *ifp1;
	FILE *ifp2;

	ifp1 = fopen("cg.lammpstrj.physical","r");
	ifp2 = fopen("dump.lammpstrj.virtrual","r");

	double x[2];
	double y[2];
	double z[2];

	double fx[2];
	double fy[2];
	double fz[2];

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

	for(t=0;t<100;t++){

		printf("%d\n",t);

		fgets(buffer,100,ifp1);
		fgets(buffer,100,ifp2);
//		fprintf(ofp,"ITEM: TI2MESTEP\n",buffer);
		
		printf("first buffer\n");

		fgets(buffer,100,ifp1);
		fgets(buffer,100,ifp2);
//		fprintf(ofp,"%d\n",ts);

		fgets(buffer,100,ifp1);
		fgets(buffer,100,ifp2);
//		fprintf(ofp,"ITEM: NUMBER OF ATOMS\n");

		fgets(buffer,100,ifp1);
		fgets(buffer,100,ifp2);
//		fprintf(ofp,"250\n");

		fgets(buffer,100,ifp1);
		fgets(buffer,100,ifp2);
//		fprintf(ofp,"ITEM: BOX BOUNDS pp pp pp\n");

		fgets(buffer,100,ifp1);
		fgets(buffer,100,ifp2);
//		fprintf(ofp,"0 31.66\n");

		fgets(buffer,100,ifp1);
		fgets(buffer,100,ifp2);
//		fprintf(ofp,"0 31.66\n");

		fgets(buffer,100,ifp1);
		fgets(buffer,100,ifp2);
//		fprintf(ofp,"0 31.66\n");

		fgets(buffer,100,ifp1);
		fgets(buffer,100,ifp2);
//		fprintf(ofp,"ITEM: ATOMS id type x y z fx fy fz \n");
//		printf("%f\n",box);

		printf("after header\n");

		for(i=0;i<250;i++){
			fscanf(ifp1,"%d %d %lf %lf %lf %lf %lf %lf\n",&id,&type,&x[0],&y[0],&z[0],&fx[0],&fy[0],&fz[0]);
			fscanf(ifp2,"%d %d %lf %lf %lf %lf %lf %lf\n",&id,&type,&x[1],&y[1],&z[1],&fx[1],&fy[1],&fz[1]);


			printf("%lf %lf %lf\n",(x[1]-x[0]),(y[1]-y[0]),(z[0]-z[1]));
		}
	}
	return 0;
}
