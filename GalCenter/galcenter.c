#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include<unistd.h>
#include<sys/types.h>
#include"Memory.h"
/*
#include "tree.h"
*/
#include "header.h"
#include "ramses.h"
#include "tree.h"
#include "defs.h"
#include "hfind.h"
#include "galcenter.h"

#define REQUIRED_OMEGA_M 0.287845f
#define REQUIRED_OMEGA_B 0.047143f
#define REQUIRED_H0 69.33f
#define REQUIRED_OMEGA_L 0.712155f

FoFTPtlStruct *rbuffer;


size_t nbuffer=3000000;



size_t ngalmax = 10000000;
size_t ngal;
GalCenter *output;


int sortcenter(const void *a, const void *b){
	GalCenter *aa = (GalCenter*)a;
	GalCenter *bb = (GalCenter*)b;
	if(aa->haloid > bb->haloid) return 1;
	else if(aa->haloid < bb->haloid) return -1;
	else {
		if(aa->subgalid > bb->subgalid) return 1;
		else if(aa->subgalid<bb->subgalid) return -1;
		else return 0;
	}
}



void FREAD(FoFTPtlStruct *p, SubInfo *galinfo, FILE *fp){
	size_t i;
	if(galinfo->npall > nbuffer) {
		rbuffer=(FoFTPtlStruct*)Realloc(rbuffer, sizeof(FoFTPtlStruct)*galinfo->npall);
		nbuffer = galinfo->npall;
	}
	DmType *dm=(DmType*)rbuffer;
	GasType *gas = (GasType*)rbuffer;
	SinkType *sink = (SinkType*)rbuffer;
	StarType *star = (StarType*)rbuffer;
	fread(dm, sizeof(DmType), galinfo->npdm,fp);
	for(i=0;i<galinfo->npdm;i++) {
		p->type = TYPE_DM;
		p->x = dm[i].x;
		p->y = dm[i].y;
		p->z = dm[i].z;
		p->vx = dm[i].vx;
		p->vy = dm[i].vy;
		p->vz = dm[i].vz;
		p->mass = dm[i].mass;
//		p->link02 =  0.2*pow(dm[i].mass/2.7755e11L/omep, 0.33333333333333333333L);
//		(p++)->p.dm = dm[i];
		p++;
	}
	fread(gas, sizeof(GasType), galinfo->npgas,fp);
	for(i=0;i<galinfo->npgas;i++) {
		p->type = TYPE_GAS;
		p->x = gas[i].x;
		p->y = gas[i].y;
		p->z = gas[i].z;
		p->vx = gas[i].vx;
		p->vy = gas[i].vy;
		p->vz = gas[i].vz;
		p->mass = gas[i].mass;
//		p->link02 =  0.2*pow(gas[i].mass/2.7755e11L/omep, 0.33333333333333333333L);
//		(p++)->p.gas = gas[i];
		p++;
	}
	fread(sink, sizeof(SinkType), galinfo->npsink,fp);
	for(i=0;i<galinfo->npsink;i++) {
		p->type = TYPE_SINK;
		p->x = sink[i].x;
		p->y = sink[i].y;
		p->z = sink[i].z;
		p->vx = sink[i].vx;
		p->vy = sink[i].vy;
		p->vz = sink[i].vz;
		p->mass = sink[i].mass;
//		p->link02 =  0.2*pow(sink[i].mass/2.7755e11L/omep, 0.33333333333333333333L);
//		(p++)->p.sink = sink[i];
		p++;
	}
	fread(star, sizeof(StarType), galinfo->npstar,fp);
	for(i=0;i<galinfo->npstar;i++) {
		p->type = TYPE_STAR;
		p->x = star[i].x;
		p->y = star[i].y;
		p->z = star[i].z;
		p->vx = star[i].vx;
		p->vy = star[i].vy;
		p->vz = star[i].vz;
		p->mass = star[i].mass;
//		p->link02 =  0.2*pow(star[i].mass/2.7755e11L/omep, 0.33333333333333333333L);
//		(p++)->p.star = star[i];
		p++;
	}
}
/*
typedef struct Halo{
	size_t np;
	POSTYPE x,y,z;
	float mass, mstar,mgas,mdm,msink;
	float vx,vy,vz;
} Halo;
*/
FoFTPtlStruct *bp,*sbp,*rbp;

int myid,nid;
long long sgalnum,rgalnum, numhalo, numgal;
lint *ptl2halonum;
int mpeak;
int ready=READY,writing=WRITING;
int snp,rnp;
FILE *wrp,*wbp, *wlist;
/*                              */
float amax,a,rng,size,hubble;
int ng,nspace;
float omep,omepb,omeplam,epsilon;
int nx,ny,nz;
/*                              */
double onesolarmass=1.989E33L;
double com2real,real2com,potentfact;
double pntmass,r1kineticfact,r2kineticfact;
float acoeff[8];

void write_data(GalCenter);

void pspline_(void);
long long iii;

int main(int argc, char *argv[]) {
	MPI_Status mstatus,cstatus;
	MPI_Request request;
	char filetag[80],filename[80];
	int i,j,snp,rnp,flag;
	int finish=-99999;
	int snend,rnend;
	int nstep;
	int src,dest;
	int ii;
	HaloInfo haloinfo;
	SubInfo galinfo;
	/*          */
	/*
	struct HaloPos *csx,*crx;
	struct HaloVel *csvx,*crvx;
	float *floatsx, *floatrx;
	int *swap;
	*/
	float tmp_tidal,tmp;
	/*          */
	FILE *fp,*rhfp,*rtidalf;
	char rheader[80],rfilename[80],wfilename[80], wlistname[190];
	char dir[190];
	GalCenter *galcenters;





	(void )MPI_Init(&argc,&argv);
	(void )MPI_Comm_size(MPI_COMM_WORLD,&nid);
	(void )MPI_Comm_rank(MPI_COMM_WORLD,&myid);


	if(argc != 2){
		if(myid == motherrank)
		fprintf(stderr,".... galfind.exe [nstep] \n");
		exit(1);
	}
	nstep = atoi(argv[1]);

	sprintf(dir,"./FoF_Data/FoF.%.5d/",nstep);

	if(Make_Total_Memory(NMEG*MEGABYTE)==0){
		fprintf(stderr,"P%d Error initializing mem %ldMB\n",myid,NMEG);
	}

	/*
	MPI_Bcast(m_tidal,NUM_MASS,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(r_tidal,NUM_MASS,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	*/
	sprintf(rheader,"%s/FoF_halo_cat.%.5d",dir,nstep);
	if(myid ==0) {
		fprintf(stderr,"opening header file: %s", rheader);
		rhfp = fopen(rheader,"r");
		{
			float bias;
			float astep,anow;
			int npower;
         	fread(&size,sizeof(float),1,rhfp);
		 	fread(&hubble,sizeof(float),1,rhfp);
			fread(&omep,sizeof(float),1,rhfp);
			fread(&omepb,sizeof(float),1,rhfp);
			fread(&omeplam,sizeof(float),1,rhfp);
			omep = REQUIRED_OMEGA_M;
			omepb = REQUIRED_OMEGA_B;
			omeplam = REQUIRED_OMEGA_L;
			hubble = REQUIRED_H0;
			fread(&amax,sizeof(float),1,rhfp);
			fread(&anow,sizeof(float),1,rhfp);
			ng = nx;
			a = anow;
		}
		fprintf(stderr," well open header file\n");
	}
	MPI_Bcast(&ng,1,MPI_INT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&nspace,1,MPI_INT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&omep,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&omeplam,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&hubble,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&size,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&amax,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&a,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	if(myid==motherrank){
	}
	rng = ng = nspace = 1; /* obsolite variables */
	nx = ny = nz = ng;
	epsilon = (nspace*0.1)*(nspace*0.1);
	if(myid==-1){
		int kkk = 1;
		while(kkk){
			kkk = 1;
		}
	}
	if(0){ 
		pid_t pid = getpid();
		char outf[9000]; 
		sprintf(outf,"CHECK/indcheck.%.5d.%.5d", nstep, myid); 
		FILE *wp = fopen(outf,"w"); 
		fprintf(wp, "pid %d \n", pid);
		fflush(wp); 
		fclose(wp);
	}

	acoeff[0] = 1.;
#ifdef DEBUG
	if(myid==0) fprintf(stderr,"Well calculating physical parameters size=%g omepb=%g potentfact=%g\n",
			size,omepb,potentfact);
#endif
	snend = rnend = INIT_NP;
	sbp = (FoFTPtlStruct*)Malloc(sizeof(FoFTPtlStruct)*snend,PPTR(sbp));
	if(myid==motherrank) {
		rbp = (FoFTPtlStruct*)Malloc(sizeof(FoFTPtlStruct)*INIT_NP,PPTR(rbp));
	}

	/*
	for(ii=4;ii<100;ii++){
	*/
	ii = 0;
	{	
		if(myid == motherrank) {
			char backfile[190];

			numhalo = sgalnum = rgalnum = 0;
			sprintf(rfilename,"%sGALFIND.DATA.%.5d",dir,nstep);
			sprintf(wfilename,"%sGALFIND.CENTER.%.5d",dir,nstep);
			fprintf(stdout,"Opening %s file for output\n",rfilename);
			FILE *frp = fopen(rfilename,"r");

			output = (GalCenter*)malloc(sizeof(GalCenter)*ngalmax);
			ngal = 0;


			rbuffer = (FoFTPtlStruct*)Malloc(sizeof(FoFTPtlStruct)*nbuffer,PPTR(rbuffer));
			while(fread(&haloinfo,sizeof(HaloInfo),1,frp) == 1){
				snp = haloinfo.npall; snend = MAX(snend,snp);
				if(snp == snend) sbp = Realloc(sbp,sizeof(FoFTPtlStruct)*snend);

				int isub;
				for(isub=0;isub<haloinfo.nsub;isub++){
					fread(&galinfo,sizeof(SubInfo), 1, frp);
					FREAD(sbp,&galinfo,frp);

					snp = galinfo.npall;
					if(numhalo%1000==0) printf("Now passing through %d:   %d\n",numhalo, snend);
					if(snp >= 5000000) printf("Now passing throught %d:   %d\n",numhalo, haloinfo.npall);

					GalCenter sgalcenter;
					sgalcenter.haloid = numhalo;
					sgalcenter.subgalid = isub;
//					if(numhalo != 5147565) goto there;
					do {
						MPI_Probe(MPI_ANY_SOURCE,READY,MPI_COMM_WORLD,&mstatus);
						src = mstatus.MPI_SOURCE;
						dest = mstatus.MPI_SOURCE;
						MPI_Recv(&ready,1,MPI_INT,src,READY,MPI_COMM_WORLD,&cstatus);
						if(ready == READY){
							sgalnum++;
							MPI_Send(&snp,1,MPI_INT,dest,NP_TAG, MPI_COMM_WORLD);
							MPI_Send(&sgalcenter,sizeof(GalCenter),MPI_BYTE,dest,NP_TAG, MPI_COMM_WORLD);
							if(snp*sizeof(FoFTPtlStruct)>1500000000L) BIG_MPI_Send(sbp,snp*sizeof(FoFTPtlStruct),MPI_BYTE,dest,R_TAG,MPI_COMM_WORLD);
							else MPI_Send(sbp,snp*sizeof(FoFTPtlStruct),MPI_BYTE,dest,R_TAG,MPI_COMM_WORLD);
						}
						else {
							rgalnum++;
							MPI_Probe(src,NP_TAG,MPI_COMM_WORLD,&cstatus);
							GalCenter rgalcenter;
							MPI_Recv(&rgalcenter,sizeof(GalCenter),MPI_BYTE,src,NP_TAG,MPI_COMM_WORLD,&cstatus);
							/*
							output[ngal++] = rgalcenter;
							if(ngal >=ngalmax) {
								ngalmax += 4000000;
								output =(GalCenter*)realloc(output, sizeof(GalCenter)*ngalmax);
							}
							*/
							write_data(rgalcenter);
						}
					} while(ready != READY); 
there:
					continue;
				}
				numhalo ++;
			}
			fclose(frp);

			printf("complete................... \n");
			fflush(stdout);
			j = 0;
			for(iii=rgalnum+1;iii<=sgalnum;){
				MPI_Probe(MPI_ANY_SOURCE,READY,MPI_COMM_WORLD,&mstatus);
				MPI_Recv(&ready,1,MPI_INT,mstatus.MPI_SOURCE,READY,
						MPI_COMM_WORLD,&cstatus);
				if(ready == WRITING){
					iii++;
					MPI_Probe(mstatus.MPI_SOURCE,NP_TAG,MPI_COMM_WORLD,&cstatus);
					GalCenter rgalcenter;
					MPI_Recv(&rgalcenter,sizeof(GalCenter),MPI_BYTE,mstatus.MPI_SOURCE,
							NP_TAG,MPI_COMM_WORLD, &cstatus);
					/*
					output[ngal++] = rgalcenter;
					if(ngal >=ngalmax) {
						ngalmax += 4000000;
						output =(GalCenter*)realloc(output, sizeof(GalCenter)*ngalmax);
					}
					*/
					write_data( rgalcenter);
				}
				else {
					j++;
					MPI_Send(&finish,1,MPI_INT,mstatus.MPI_SOURCE, NP_TAG,MPI_COMM_WORLD);
				}
			}
			for(i=1;i<nid-j;i++){
				MPI_Recv(&ready,1,MPI_INT,MPI_ANY_SOURCE,READY, MPI_COMM_WORLD,&mstatus);
				MPI_Send(&finish,1,MPI_INT,mstatus.MPI_SOURCE, NP_TAG,MPI_COMM_WORLD);
			}

			qsort(output, ngal, sizeof(GalCenter), sortcenter);
			wrp = fopen(wfilename,"w");
			fwrite(output, sizeof(GalCenter), ngal, wrp);
			fclose(wrp);
		} /* end of "myid" for */
		else {
			ready = READY;
			MPI_Send(&ready,1,MPI_INT,motherrank,READY,MPI_COMM_WORLD);
			MPI_Recv(&snp,1,MPI_INT,motherrank,NP_TAG,MPI_COMM_WORLD,&cstatus);
			while(snp != finish){
				GalCenter rgalcenter,tmpgal;
				long numstack;
				MPI_Recv(&rgalcenter,sizeof(GalCenter),MPI_BYTE,motherrank,NP_TAG,MPI_COMM_WORLD,&cstatus);
				snend = MAX(snend,snp);
				sbp = Realloc(sbp,sizeof(FoFTPtlStruct)*snend);
				if(snp*sizeof(FoFTPtlStruct)>1500000000L)  BIG_MPI_Recv(sbp,sizeof(FoFTPtlStruct)*snp,MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
				else MPI_Recv(sbp,sizeof(FoFTPtlStruct)*snp,MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
	
				if(0){
					char outf[9000];
					sprintf(outf,"CHECK/indcheck.%.5d.%.5d", nstep, myid);
					FILE *wp = fopen(outf,"a");
					fprintf(wp, "%ld %ld %d", rgalcenter.haloid,rgalcenter.subgalid, snp);
					fflush(wp);
					fclose(wp);
				}
				numstack = CurMemStack();
				find_centers(&rgalcenter, sbp,snp);
				InitialOldMemStack(numstack);
				/*
				printf("End.initial send"); fflush(stdout);
				*/
	
				ready = WRITING;
				MPI_Issend(&ready,1,MPI_INT,motherrank,READY,MPI_COMM_WORLD, &request);
				MPI_Wait(&request,&mstatus);
				MPI_Send(&rgalcenter,sizeof(GalCenter),MPI_BYTE,motherrank,NP_TAG, MPI_COMM_WORLD);
				/*
				printf("sent\n "); fflush(stdout);
				*/
				if(0){
					char outf[9000];
					sprintf(outf,"CHECK/indcheck.%.5d.%.5d", nstep, myid);
					FILE *wp = fopen(outf,"a");
					fprintf(wp, "   returned: %d\n",snp);
					fflush(wp);
					fclose(wp);
				}
	
				ready = READY;
				MPI_Send(&ready,1,MPI_INT,motherrank,READY,MPI_COMM_WORLD);
				MPI_Recv(&snp,1,MPI_INT,motherrank,NP_TAG,MPI_COMM_WORLD,&cstatus);
				if(snend > 2*INIT_NP) snend = 2*INIT_NP;
			}
			if(0){
				char outf[9000];
				sprintf(outf,"CHECK/indcheck.%.5d.%.5d", nstep, myid);
				FILE *wp = fopen(outf,"a");
				fprintf(wp, "End of Job\n");
					fflush(wp);
				fclose(wp);
			}
		}
	}
	if(myid == motherrank) fprintf(stdout,"End of calculation. Closing\n");
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
void bin2ascii(struct PtlPos *,struct PtlVel *,lint );
void write_data(GalCenter rgalcenter){
	output[ngal++] = rgalcenter;
	if(ngal >=ngalmax) {
		ngalmax += 4000000;
		output =(GalCenter*)realloc(output, sizeof(GalCenter)*ngalmax);
	}
}

