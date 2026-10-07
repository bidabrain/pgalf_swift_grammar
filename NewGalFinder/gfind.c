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
static int tidal_start=1;
float m_tidal[NUM_MASS],r_tidal[NUM_MASS];


FoFTPtlStruct *rbuffer;
size_t nbuffer=3000000;

#define REQUIRED_OMEGA_M 0.287845f
#define REQUIRED_OMEGA_B 0.047143f
#define REQUIRED_H0 69.33f
#define REQUIRED_OMEGA_L 0.712155f

void FREAD(FoFTPtlStruct *p, HaloQ *haloq, FILE *fp, float omep){
	const float omega_safe = (omep > 0.f) ? omep : REQUIRED_OMEGA_M;
	size_t i;
	if(haloq->np > nbuffer) {
		rbuffer=(FoFTPtlStruct*)Realloc(rbuffer, sizeof(FoFTPtlStruct)*haloq->np);
		nbuffer = haloq->np;
	}
	DmType *dm=(DmType*)rbuffer;
	GasType *gas = (GasType*)rbuffer;
	SinkType *sink = (SinkType*)rbuffer;
	StarType *star = (StarType*)rbuffer;
	fread(dm, sizeof(DmType), haloq->npdm,fp);
#ifdef DEBUG
//	printf("reading data: %ld %ld %ld %ld \n", haloq->npdm, haloq->npgas, haloq->npsink,
//			haloq->npstar);
#endif
	for(i=0;i<haloq->npdm;i++) {
		p->type = TYPE_DM;
		p->x = dm[i].x;
		p->y = dm[i].y;
		p->z = dm[i].z;
		p->vx = dm[i].vx;
		p->vy = dm[i].vy;
		p->vz = dm[i].vz;
		p->mass = dm[i].mass;
		p->link02 =  0.2*pow(dm[i].mass/2.7755e11L/omega_safe, 0.33333333333333333333L);
		(p++)->p.dm = dm[i];
	}
	fread(gas, sizeof(GasType), haloq->npgas,fp);
	for(i=0;i<haloq->npgas;i++) {
		p->type = TYPE_GAS;
		p->x = gas[i].x;
		p->y = gas[i].y;
		p->z = gas[i].z;
		p->vx = gas[i].vx;
		p->vy = gas[i].vy;
		p->vz = gas[i].vz;
		p->mass = gas[i].mass;
		p->link02 =  0.2*pow(gas[i].mass/2.7755e11L/omega_safe, 0.33333333333333333333L);
		(p++)->p.gas = gas[i];
	}
	fread(sink, sizeof(SinkType), haloq->npsink,fp);
	for(i=0;i<haloq->npsink;i++) {
		p->type = TYPE_SINK;
		p->x = sink[i].x;
		p->y = sink[i].y;
		p->z = sink[i].z;
		p->vx = sink[i].vx;
		p->vy = sink[i].vy;
		p->vz = sink[i].vz;
		p->mass = sink[i].mass;
		p->link02 =  0.2*pow(sink[i].mass/2.7755e11L/omega_safe, 0.33333333333333333333L);
		(p++)->p.sink = sink[i];
	}
	fread(star, sizeof(StarType), haloq->npstar,fp);
	for(i=0;i<haloq->npstar;i++) {
		p->type = TYPE_STAR;
		p->x = star[i].x;
		p->y = star[i].y;
		p->z = star[i].z;
		p->vx = star[i].vx;
		p->vy = star[i].vy;
		p->vz = star[i].vz;
		p->mass = star[i].mass;
		p->link02 =  0.2*pow(star[i].mass/2.7755e11L/omega_safe, 0.33333333333333333333L);
		(p++)->p.star = star[i];
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
long long shalonum,rhalonum, numhalo;
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

void write_data(FoFTPtlStruct *,lint*,int);

int physical_parameters(void);
int subhalo_den(FoFTPtlStruct *, lint,lint *);
void pspline_(void);
void i_potent_spline(void);
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
	HaloQ haloq;
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





	MPI_Init(&argc,&argv);

	MPI_Comm_size(MPI_COMM_WORLD,&nid);

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);


#ifdef NSPLIT
	if(argc != 3){
		if(myid == motherrank)
		fprintf(stderr,".... hfind.exe [nstep] [mid]\n");
		exit(1);
	}
#else
	if(argc != 2 && argc != 4){
		if(myid == motherrank)
		fprintf(stderr,".... hfind.exe [nstep] \n");
		exit(1);
	}
#endif
	nstep = atoi(argv[1]);
	long ioffset, joffset;
	if(argc ==4){
		ioffset = atol(argv[2]);
		joffset = atol(argv[3]);
	}
	

	sprintf(dir,"./FoF_Data/FoF.%.5d/",nstep);

#ifdef DEBUG
	printf("Data directory: %s\n", dir);
#endif

	INT8 msize = NMEG * MEGABYTE;
	if(myid ==0) msize = 80000*MEGABYTE;
	if(Make_Total_Memory(msize )==0){
		fprintf(stderr,"P%d Error initializing mem %ldMB\n",myid,NMEG);
	}

	i_potent_spline();
	if(myid == motherrank && tidal_start == 1){
		/*
		fprintf(stderr,"reading tidal.dat");
		if((rtidalf = fopen("tidal.dat","r"))==NULL){
			fprintf(stderr,"Error opening tidal.dat\n");
			exit(9);
		}
			
		for(i=0;i<NUM_MASS;i++){
			fscanf(rtidalf,"%g %g %g\n",&m_tidal[i],&r_tidal[i],&tmp);
		}
		fclose(rtidalf);
		tidal_start = 0;
		fprintf(stderr,"read tidal.dat");
		*/
	}
	/*
	MPI_Bcast(m_tidal,NUM_MASS,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(r_tidal,NUM_MASS,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	*/
#ifdef NSPLIT
	ii = 0;
	sprintf(rheader,"%s/FoF_halo_cat.%.5d.%.5d",dir,nstep,ii);
#else
	sprintf(rheader,"%s/FoF_halo_cat.%.5d",dir,nstep);
#endif
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
	(void) physical_parameters();
#ifdef DEBUG
	if(myid==0) fprintf(stderr,"Well calculating physical parameters size=%g omepb=%g potentfact=%g omep= %g\n",
			size,omepb,potentfact, omep);
#endif
#ifdef PMFORCE
    pspline_();
#endif
	snend = rnend = INIT_NP;
	sbp = (FoFTPtlStruct*)Malloc(sizeof(FoFTPtlStruct)*snend,PPTR(sbp));
	ptl2halonum = (lint *)Malloc(sizeof(lint)*INIT_NP,PPTR(ptl2halonum));
	if(myid==motherrank) {
		rbp = (FoFTPtlStruct*)Malloc(sizeof(FoFTPtlStruct)*INIT_NP,PPTR(rbp));
	}

	/*
	for(ii=4;ii<100;ii++){
	*/
#ifdef NSPLIT
	ii = atoi(argv[2]);
#else
	ii = 0;
#endif
	{	
		if(myid == motherrank) {
			char backfile[190];
#ifdef NSPLIT
			sprintf(backfile,"%sbackground_ptl.%.5d.%.5d",dir,nstep,ii);
			wbp = fopen(backfile,"w");

			numhalo = shalonum = rhalonum = 0;
			sprintf(wfilename,"%sGALFIND.DATA.%.5d.%.5d",dir,nstep,ii);
			sprintf(wlistname,"%sGALCATALOG.LIST.%.5d.%.5d",dir,nstep,ii);
			fprintf(stdout,"Opening %s file for output\n",wfilename);
			wrp = fopen(wfilename,"w");
			wlist = fopen(wlistname,"w");

			sprintf(rfilename,"%sFoF_member_particle.%.5d.%.5d",dir,nstep,ii);
			fp = fopen(rfilename,"r");
			if(ii!=0) {
				sprintf(rheader,"%sFoF_halo_cat.%.5d.%.5d",dir,nstep,ii);
				rhfp = fopen(rheader,"r");
			}
#else
			sprintf(backfile,"%sbackground_ptl.%.5d",dir,nstep);
			wbp = fopen(backfile,"w");

			numhalo = shalonum = rhalonum = 0;
			sprintf(wfilename,"%sGALFIND.DATA.%.5d",dir,nstep);
			sprintf(wlistname,"%sGALCATALOG.LIST.%.5d",dir,nstep);
			fprintf(stdout,"Opening %s file for output\n",wfilename);
			wrp = fopen(wfilename,"w");
			wlist = fopen(wlistname,"w");

			sprintf(rfilename,"%sFoF_member_particle.%.5d",dir,nstep);
			fp = fopen(rfilename,"r");
#endif
			rbuffer = (FoFTPtlStruct*)Malloc(sizeof(FoFTPtlStruct)*nbuffer,PPTR(rbuffer));
			if(0){
				HaloQ rbuff[10000];
				long offset = 0;
				for(i=0;i<190;i++) {
					int np;
					np = fread(rbuff, sizeof(HaloQ), 10000,rhfp);
					for(j=0;j<np;j++){
						offset += sizeof(DmType)*rbuff[j].npdm;
						offset += sizeof(GasType)*rbuff[j].npgas;
						offset += sizeof(SinkType)*rbuff[j].npsink;
						offset += sizeof(StarType)*rbuff[j].npstar;
					}
					numhalo += np;
				}
				if(0){
					int np = fread(rbuff, sizeof(HaloQ), 0, rhfp);
					for(j=0;j<np;j++){
						offset += sizeof(DmType)*rbuff[j].npdm;
						offset += sizeof(GasType)*rbuff[j].npgas;
						offset += sizeof(SinkType)*rbuff[j].npsink;
						offset += sizeof(StarType)*rbuff[j].npstar;
					}
					numhalo += np;
				}

				printf("now jumping to %ld\n",offset);
				fseek(fp, offset, SEEK_CUR);

			}
			if(argc==4){
				fseek(rhfp,ioffset, SEEK_SET);
				fseek(fp,joffset, SEEK_SET);
			}
			while(fread(&haloq,sizeof(HaloQ),1,rhfp) == 1){
				size_t qoffset = ftell(rhfp) - sizeof(HaloQ);
				snp = haloq.np; snend = MAX(snend,snp);
				if(snp == snend) sbp = Realloc(sbp,sizeof(FoFTPtlStruct)*snend);
				size_t soffset = ftell(fp);
{

}
				FREAD(sbp,&haloq,fp, omep);
				numhalo ++;
				if(numhalo%1000==0) printf("Now passing through %lld:   %d for snp= %d offsets= %ldL %ldL\n",numhalo, snend, snp, qoffset, soffset);
				/*
				if(snp > 5 && snp > 100000 && numhalo > 1278000) 
				if(snp > 5 && snp > 1000000) ===> fails
				if(snp > 5 && snp <= 1000000 && snp > 100000) => fails at numhalo>1423000
				if(snp <= 1000000 && snp > 100000 && numhalo >= 1423747)
				if(snp == 5897207 )
				if(snp == 1209448 )
				if(snp >= 30 && snp < 100000 && numhalo > 6000000)
				*/
				/*
				*/
				
//				if(snp == 165557336)
//				if(snp >= 10000 )
				if(snp >= 30 )
				
				{
//					if(snp >= 100000) 
						LOGPRINT("Now passing through %lld:   %zu with offsets= %ldL %ldL\n",
								numhalo, haloq.np, qoffset, soffset);
					do {
						MPI_Probe(MPI_ANY_SOURCE,READY,MPI_COMM_WORLD,&mstatus);
						src = mstatus.MPI_SOURCE;
						dest = mstatus.MPI_SOURCE;
						MPI_Recv(&ready,1,MPI_INT,src,READY,MPI_COMM_WORLD,&cstatus);
						if(ready == READY){
							shalonum++;
							MPI_Send(&snp,1,MPI_INT,dest,NP_TAG, MPI_COMM_WORLD);
							MPI_Send(&numhalo,1,MPI_LONG_LONG,dest,NP_TAG, MPI_COMM_WORLD);
							if(snp*sizeof(FoFTPtlStruct)>1500000000L) BIG_MPI_Send(sbp,snp*sizeof(FoFTPtlStruct),MPI_BYTE,dest,R_TAG,MPI_COMM_WORLD);
							else MPI_Send(sbp,snp*sizeof(FoFTPtlStruct),MPI_BYTE,dest,R_TAG,MPI_COMM_WORLD);
#ifdef DEBUG
//							fprintf(stderr,"P0 sending to %d with %d particles\n",mstatus.MPI_SOURCE,snp); 
//							fflush(stderr);
#endif
						}
						else {
							rhalonum++;
							MPI_Probe(src,NP_TAG,MPI_COMM_WORLD,&cstatus);
							MPI_Recv(&rnp,1,MPI_INT,src,NP_TAG,MPI_COMM_WORLD,&cstatus);
							/*
							fprintf(stderr,"receiving from %d with %d particles\n",
								src,rnp); fflush(stderr);
								*/
							if(rnp > rnend) {
								rnend = MAX(rnend,rnp);
								rbp = Realloc(rbp,sizeof(FoFTPtlStruct)*rnend);
								ptl2halonum = Realloc(ptl2halonum,sizeof(lint)*rnend);
							}
							MPI_Probe(src,R_TAG,MPI_COMM_WORLD,&cstatus);
							if(rnp*sizeof(FoFTPtlStruct)>1500000000L) BIG_MPI_Recv(rbp,sizeof(FoFTPtlStruct)*rnp,MPI_BYTE,src,R_TAG,MPI_COMM_WORLD,&cstatus);
							else  MPI_Recv(rbp,sizeof(FoFTPtlStruct)*rnp,MPI_BYTE,src,R_TAG,MPI_COMM_WORLD,&cstatus);
							MPI_Probe(src,PT2H_TAG,MPI_COMM_WORLD,&cstatus);
							MPI_Recv(ptl2halonum,rnp,MPI_INT,src, PT2H_TAG, MPI_COMM_WORLD,&cstatus);
							MPI_Probe(src,MPEAK_TAG,MPI_COMM_WORLD,&cstatus);
							MPI_Recv(&mpeak,1,MPI_INT,src,MPEAK_TAG, MPI_COMM_WORLD,&cstatus);
							if(mpeak < 0 || mpeak > rnp) mpeak = 0;
							write_data(rbp,ptl2halonum,rnp);
#ifdef DEBUG
//							printf("P0 %d received with mpeak=%d\n",rnp,mpeak);
#endif
						}
					} while(ready != READY); 
					/*
					if(rnend > 2*INIT_NP) rnend = 2*INIT_NP;
					if(snend > 2*INIT_NP) snend = 2*INIT_NP;
					*/
				}
				if(argc==4) break;
			}
			fclose(fp);
			fclose(rhfp);

			printf("complete................... \n");
			fflush(stdout);
			j = 0;
			for(iii=rhalonum+1;iii<=shalonum;){
				MPI_Probe(MPI_ANY_SOURCE,READY,MPI_COMM_WORLD,&mstatus);
				MPI_Recv(&ready,1,MPI_INT,mstatus.MPI_SOURCE,READY,
						MPI_COMM_WORLD,&cstatus);
				if(ready == WRITING){
					iii++;
					MPI_Probe(mstatus.MPI_SOURCE,NP_TAG,MPI_COMM_WORLD,&cstatus);
					MPI_Recv(&rnp,1,MPI_INT,mstatus.MPI_SOURCE,NP_TAG,MPI_COMM_WORLD, &cstatus);
					rnend = MAX(rnend,rnp);
					rbp = Realloc(rbp,sizeof(FoFTPtlStruct)*rnend);
					ptl2halonum = Realloc(ptl2halonum,sizeof(lint)*rnend);
					MPI_Probe(mstatus.MPI_SOURCE,R_TAG,MPI_COMM_WORLD, &cstatus);
					if(rnp*sizeof(FoFTPtlStruct)>1500000000L) BIG_MPI_Recv(rbp,sizeof(FoFTPtlStruct)*rnp,MPI_BYTE,mstatus.MPI_SOURCE,R_TAG,MPI_COMM_WORLD,&cstatus);
					else MPI_Recv(rbp,sizeof(FoFTPtlStruct)*rnp,MPI_BYTE,mstatus.MPI_SOURCE,R_TAG,MPI_COMM_WORLD,&cstatus);
					MPI_Probe(mstatus.MPI_SOURCE,PT2H_TAG,MPI_COMM_WORLD, &cstatus);
					MPI_Recv(ptl2halonum,rnp,MPI_INT,mstatus.MPI_SOURCE,PT2H_TAG, MPI_COMM_WORLD,&cstatus);
					MPI_Probe(mstatus.MPI_SOURCE,MPEAK_TAG,MPI_COMM_WORLD, &cstatus);
					MPI_Recv(&mpeak,1,MPI_INT,mstatus.MPI_SOURCE,MPEAK_TAG, MPI_COMM_WORLD,&cstatus);
					if(mpeak < 0 || mpeak > rnp) mpeak = 0;
					if(1) fprintf(stdout,"receiving from %d with %d particles : %lld : %lld\n",
								mstatus.MPI_SOURCE,rnp,iii,shalonum); fflush(stdout);
					write_data(rbp,ptl2halonum,rnp);
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
			fclose(wrp);
			fclose(wlist);
			fclose(wbp);
		} /* end of "myid" for */
		else {
			ready = READY;
			MPI_Send(&ready,1,MPI_INT,motherrank,READY,MPI_COMM_WORLD);
			MPI_Recv(&snp,1,MPI_INT,motherrank,NP_TAG,MPI_COMM_WORLD,&cstatus);
			while(snp != finish){
				long long inumhalo;
				long numstack;
				MPI_Recv(&inumhalo,1,MPI_LONG_LONG,motherrank,NP_TAG,MPI_COMM_WORLD,&cstatus);
				snend = MAX(snend,snp);
#ifdef DEBUG
//				printf("P%d: receiving particles %d ",myid,snp); fflush(stdout);
#endif
				sbp = Realloc(sbp,sizeof(FoFTPtlStruct)*snend);
				ptl2halonum = Realloc(ptl2halonum,sizeof(lint)*snend);
				if(snp*sizeof(FoFTPtlStruct)>1500000000L)  BIG_MPI_Recv(sbp,sizeof(FoFTPtlStruct)*snp,MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
				else MPI_Recv(sbp,sizeof(FoFTPtlStruct)*snp,MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
	
				if(0){
					char outf[9000];
					sprintf(outf,"CHECK/indcheck.%.5d.%.5d", nstep, myid);
					FILE *wp = fopen(outf,"a");
					fprintf(wp, "%lld %d", inumhalo, snp);
					fflush(wp);
					fclose(wp);
				}
				numstack = CurMemStack();
				mpeak = subhalo_den(sbp,snp,ptl2halonum);
				InitialOldMemStack(numstack);
#ifdef DEBUG
//				printf("P%d: End.initial send", myid); fflush(stdout);
#endif
	
				ready = WRITING;
				MPI_Issend(&ready,1,MPI_INT,motherrank,READY,MPI_COMM_WORLD, &request);
				MPI_Wait(&request,&mstatus);
				MPI_Send(&snp,1,MPI_INT,motherrank,NP_TAG,MPI_COMM_WORLD);
				if(snp*sizeof(FoFTPtlStruct)>1500000000L) BIG_MPI_Send(sbp,sizeof(FoFTPtlStruct)*snp,MPI_BYTE,motherrank,R_TAG, MPI_COMM_WORLD);
				else MPI_Send(sbp,sizeof(FoFTPtlStruct)*snp,MPI_BYTE,motherrank,R_TAG, MPI_COMM_WORLD);
				MPI_Send(ptl2halonum,snp,MPI_INT,motherrank,PT2H_TAG, MPI_COMM_WORLD);
				MPI_Send(&mpeak,1,MPI_INT,motherrank,MPEAK_TAG,MPI_COMM_WORLD);
#ifdef DEBUG
//				printf("P%d: sent\n ",myid); fflush(stdout);
#endif
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
HaloInfo GetHaloInfo(FoFTPtlStruct *bp, int np){
	int i,j,k;
	HaloInfo res;
	res.nsub = mpeak;
	res.totm = res.mdm = res.mgas = res. msink = res.mstar = 0;
	res.x = res.y = res.z = res. vx = res.vy = res.vz = 0;
	res.nstar = res.nsink = res.ngas = res.ndm = 0;
	for(i=0;i<np;i++){
		dptype mass,x,y,z,vx,vy,vz;
		if(bp[i].type == TYPE_DM) {
			res.ndm ++;
			mass = bp[i].p.dm.mass;
			res.mdm += mass;
			vx = bp[i].p.dm.vx;
			vy = bp[i].p.dm.vy;
			vz = bp[i].p.dm.vz;
		}
		else if(bp[i].type == TYPE_SINK){
			res.nsink ++;
			mass = bp[i].p.sink.mass;
			res.msink += mass;
			vx = bp[i].p.sink.vx;
			vy = bp[i].p.sink.vy;
			vz = bp[i].p.sink.vz;
		}
		else if(bp[i].type == TYPE_STAR){
			res.nstar ++;
			mass = bp[i].p.star.mass;
			res.mstar += mass;
			vx = bp[i].p.star.vx;
			vy = bp[i].p.star.vy;
			vz = bp[i].p.star.vz;
		}
		else if(bp[i].type == TYPE_GAS){
			res.ngas ++;
			mass = bp[i].p.gas.mass;
			res.mgas += mass;
			vx = bp[i].p.gas.vx;
			vy = bp[i].p.gas.vy;
			vz = bp[i].p.gas.vz;
		}
		x = bp[i].x;
		y = bp[i].y;
		z = bp[i].z;
		res.totm += mass;
		res.x += mass*x;
		res.y += mass*y;
		res.z += mass*z;
		res.vx += mass*vx;
		res.vy += mass*vy;
		res.vz += mass*vz;
	}
	res.npall  = np;
	res.x = res.x/res.totm;
	res.y = res.y/res.totm;
	res.z = res.z/res.totm;
	res.vx = res.vx/res.totm;
	res.vy = res.vy/res.totm;
	res.vz = res.vz/res.totm;
	return res;
}
void GetSubHaloInfo(FoFTPtlStruct *bp, int np, lint *rptl2halonum, int id, SubInfo *shalo, size_t *asize, void *subdata){
	int i,j,k;
	void *res;
	DmType *dm, *rdm;
	StarType *star, *rstar;
	SinkType *sink, *rsink;
	GasType *gas, *rgas;
	rdm = dm = (DmType*)Malloc(sizeof(DmType)*np,PPTR(rdm));
	rstar = star = (StarType*)Malloc(sizeof(StarType)*np,PPTR(rstar));
	rsink = sink = (SinkType*)Malloc(sizeof(SinkType)*np,PPTR(rsink));
	rgas = gas = (GasType*)Malloc(sizeof(GasType)*np,PPTR(rgas));
	shalo->npall = shalo->npdm = shalo->npstar = shalo->npsink = shalo->npgas = 0;
	shalo->totm = shalo->mdm = shalo->mstar = shalo->msink = shalo->mgas = 0;
	shalo->x = shalo->y = shalo->z = 0;
	shalo->vx = shalo->vy = shalo->vz = 0;
	for(i=0;i<np;i++){
		if(rptl2halonum[i] == id){
			dptype mass;
			if(bp[i].type == TYPE_DM) {
				*(dm++) = bp[i].p.dm;
				shalo->npdm ++;
				mass = bp[i].p.dm.mass;
				shalo->mdm += mass;
			}
			else if(bp[i].type == TYPE_SINK) {
				*(sink++) = bp[i].p.sink;
				shalo->npsink ++;
				mass = bp[i].p.sink.mass;
				shalo->msink += mass;
			}
			else if(bp[i].type == TYPE_STAR) {
				*(star++) = bp[i].p.star;
				shalo->npstar ++;
				mass = bp[i].p.star.mass;
				shalo->mstar += mass;
			}
			else if(bp[i].type == TYPE_GAS) {
				*(gas++) = bp[i].p.gas;
				shalo->npgas ++;
				mass = bp[i].p.gas.mass;
				shalo->mgas += mass;
			}
			shalo->x += mass*bp[i].x;
			shalo->y += mass*bp[i].y;
			shalo->z += mass*bp[i].z;
			shalo->vx += mass*bp[i].vx;
			shalo->vy += mass*bp[i].vy;
			shalo->vz += mass*bp[i].vz;
			shalo->totm += mass;
			shalo->npall ++;
		}
	}
	shalo->x /= shalo->totm;
	shalo->y /= shalo->totm;
	shalo->z /= shalo->totm;
	shalo->vx /= shalo->totm;
	shalo->vy /= shalo->totm;
	shalo->vz /= shalo->totm;
	size_t size = 0;
	memcpy((char*)subdata+size, rdm,  (dm-rdm)*sizeof(DmType)); size += (dm-rdm)*sizeof(DmType);
	memcpy((char*)subdata+size, rgas, (gas-rgas)*sizeof(GasType)); size += (gas-rgas)*sizeof(GasType);
	memcpy((char*)subdata+size, rsink, (sink-rsink)*sizeof(SinkType)); size += (sink-rsink)*sizeof(SinkType);
	memcpy((char*)subdata+size, rstar, (star-rstar)*sizeof(StarType)); size += (star-rstar)*sizeof(StarType);

	*asize = size;

	Free(rgas);
	Free(rsink);
	Free(rstar);
	Free(rdm);

}
void write_data(FoFTPtlStruct *ssbp,lint *rptl2halonum,int rrnp){
	FoFTPtlStruct *wr;
	int wnp;
	int i,j,k;

	if(rrnp <= 0 || ssbp == NULL || rptl2halonum == NULL) return;
	if(wrp == NULL || wbp == NULL || wlist == NULL) return;
	if(mpeak < 0 || mpeak > rrnp) mpeak = 0;

	if (mpeak > 0)
	{
		HaloInfo hinfo;

		hinfo = GetHaloInfo(ssbp,rrnp);
		hinfo.nsub = mpeak;
        fwrite(&hinfo,sizeof(HaloInfo),1,wrp);
        fwrite(&hinfo,sizeof(HaloInfo),1,wlist);
		void *subdata = (void*)Malloc(sizeof(FoFTPtlStruct)*rrnp, PPTR(subdata));

		for(i=0;i<mpeak;i++){
			SubInfo subinfo;
			size_t size;
			GetSubHaloInfo(ssbp,rrnp, rptl2halonum, i, &subinfo, &size, subdata);
			fwrite(&subinfo, sizeof(SubInfo), 1, wrp);
			fwrite(&subinfo, sizeof(SubInfo), 1, wlist);
			fwrite(subdata, sizeof(char), size, wrp);

		}
#ifndef NOBACKGROUND
		i = NOT_HALO_MEMBER;
		{
			SubInfo subinfo;
			size_t size;
			GetSubHaloInfo(ssbp,rrnp, rptl2halonum, i, &subinfo, &size, subdata);
			fwrite(&subinfo, sizeof(SubInfo), 1, wbp);
			fwrite(subdata, sizeof(char), size, wbp);
		}
#endif
		Free(subdata);
	}
#ifdef CHECK_HALO
    bin2ascii(rrx,rrvx,rrnp);
#endif
	fflush(wrp);
	fflush(wbp);

}
void minmax(struct PtlPos *,int ,float *,float *, float *,float *,float *, float *,float *);
void bin2ascii(struct PtlPos *r,struct PtlVel *vr,lint np)
{
	int i,j;
	FILE *fp,*fpt,*fptz;
	char char1;
	float minx,minz,maxx,maxz,miny,maxy,width;
	int uline,dline;
	char *ctype[] = {"yellow", "blue", "red","green","6",
	"magenta","6","7","8","5","4","3","0","yellow","blue","red",
	"magenta","6","7","8","5","4","3","0","yellow","blue","red",
	"magenta","6","7","8","5","4","3","0","yellow","blue","red",
	"magenta","6","7","8","5","4","3","0","yellow","blue","red",
	"magenta","6","7","8","5","4","3","0","yellow","blue","red",
	"magenta","6","7","8","5","4","3","0","yellow","blue","red",
	"magenta","6","7","8","5","4","3","0","yellow","blue","red" };
	fp = fopen("halo0.dat","w");
	fpt = fopen("r.sm","w");
	fptz = fopen("z.sm","w");
	fprintf(fpt,"erase \n");
	fprintf(fptz,"erase \n");
	fprintf(fpt,"ticksize 0.1 1 0.1 1 \n");
	fprintf(fptz,"ticksize 0.1 1 0.1 1 \n");
	fprintf(fpt,"data halo0.dat \n");
	fprintf(fptz,"data halo0.dat \n");
	fprintf(fpt," ptype 0 0 \n");
	fprintf(fptz," ptype 0 0 \n");
	minmax(r,np,&minx,&maxx,&miny,&maxy,&minz,&maxz,&width);
	minx = floor(minx);minz = floor(minz);width = ceil(width)+1;
	fprintf(fpt," limits %f %f %f %f \n",minx,minx+width,minz,minz+width);
	fprintf(fptz," limits %f %f %f %f \n",minx,minx+width,miny,miny+width);
	fprintf(fpt," ctype default \n");
	fprintf(fptz," ctype default \n");
	fprintf(fpt," expand 1.2 \n");
	fprintf(fptz," expand 1.2 \n");
	fprintf(fpt," box \n");
	fprintf(fptz," box \n");
	fprintf(fpt," expand 0.3 \n");
	fprintf(fptz," expand 0.3 \n");
	uline = 1;dline = 0;
	for(j=0;j<mpeak;j++){
		for(i=0;i<np;i++){
			if(ptl2halonum[i] == j){
			fprintf(fp,"%f %f %f %f %f %f \n",r[i].x,r[i].y,r[i].z,
					vr[i].vx,vr[i].vy,vr[i].vz);
			dline++;
			}
		}
		fflush(fp);
		if(dline > uline ){
		fprintf(fpt," lines %d %d \n",uline,dline);
		fprintf(fptz," lines %d %d \n",uline,dline);
		fprintf(fpt," read { x 1 y 2 z 3} \n");
		fprintf(fptz," read { x 1 y 2 z 3} \n");
  		fprintf(fpt," ctype %s \n",*(ctype+j));
  		fprintf(fptz," ctype %s \n",*(ctype+j));
		fprintf(fpt," points x z \n");
		fprintf(fptz," points x y \n");
		fflush(fpt);
		}
		uline = dline + 1;
	}
	/* ?̹? ??�� ???׶????? ??ƼŬ?̴?. */
	uline = dline;
	for(i=0;i<np;i++){
		if(ptl2halonum[i]<0){
			fprintf(fp,"%f %f %f %f %f %f \n",r[i].x,r[i].y,r[i].z,
					vr[i].vx,vr[i].vy,vr[i].vz);
			dline++;
		}
	}
	fprintf(fpt," lines %d %d \n",uline,dline);
	fprintf(fptz," lines %d %d \n",uline,dline);
	fprintf(fpt," read { x 1 y 2 z 3} \n");
	fprintf(fptz," read { x 1 y 2 z 3} \n");
  	fprintf(fpt," ctype default \n");
  	fprintf(fptz," ctype default \n");
	fprintf(fpt," points x z \n");
	fprintf(fptz," points x y \n");
	fflush(fpt);
	fflush(fptz);
	fclose(fpt);
	fclose(fptz);
	fflush(fp);
	fclose(fp);
	printf("input arbitrary value : ");
	char1 = getc(stdin);
	fflush(stdin);
}

void minmax(struct PtlPos *r,int np,float *minx,float *maxx,
		float *miny,float *maxy,float *minz,
		float *maxz,float *width){
	int i;
	*minx = *miny = *minz = 1000;
	*maxx = *maxy = *maxz = 0;
	for(i=0;i<np;i++){
		*minx = MIN(*minx,r[i].x);
		*miny = MIN(*miny,r[i].y);
		*minz = MIN(*minz,r[i].z);
		*maxx = MAX(*maxx,r[i].x);
		*maxy = MAX(*maxy,r[i].y);
		*maxz = MAX(*maxz,r[i].z);
	}
	*width = MAX(*maxx-*minx,*maxz-*minz);
	*width = MAX(*width,*maxy-*miny);
	*width = MAX(*width,4.);
}

int physical_parameters(void){
	float omepk;
	float Hsub,h0= hubble/100.;
	omepk = 1.-omep-omeplam;
	Hsub = sqrt(omep*pow3(amax/a)+omeplam+omepk*pow2(amax/a));
	/* ???????? ??�� ???? ?Ѵ?. */
	/*
	com2real = size/hubble/rng/amax*a*1.E6L*pc;
	real2com = 1.L/com2real;
  	pntmass = 3.L/8.L/pi/G*pow2(100.L*100000.L*hubble)*1.E6L*pc;
	pntmass = pntmass*pow3(size/hubble)/pow3(rng/nspace)*omep;
	potentfact = G/acoeff[0]*pntmass/com2real;
	*/
	com2real = a/h0/amax*1.E6L*pc;
	real2com = 1.L/com2real;
  	pntmass = Msun/h0;
	potentfact = G/acoeff[0]*pntmass/com2real;
	/* V(r) = \dot a/ a_max *x + a/a_max * \dot x                     */
	/*      = H /(1+z) * x + u/(1+z) * \dot a                         */
	/*      = H /(1+z) * x + u/(1+z) * H * a                          */
	/*  ?׷??? x ?? ???? comoving scale ?????? distance ?̱? ??????
	 *  x = x(pm) * L{size} / rng ( rng = one dimentional mesh size)
	 *  ?̰?, 
	 *  u = u(pm) * L{size} / rng ?̾??? ?ϴµ?, ?ùķ??̼? fda ????
	 *   ?̹?, u(pm) ?ȿ? rng ?? ??????�� ?ִ?... ?ڵ忡?? pfact ??
	 *   ???? ?????? astep ?? ?־??? ?ϴµ?, rng ?? ????��?ִ?.
	 *   ?̰?�� velocity ?? ?̹? rng ?? ??????�� ?ִٴ? ??�� ???Ѵ?. ^^*/
	/* 10.24*h^-1/(amax/a)*100km/s/Mpc*h/pixel_size*(1+z)^1.5 
	r1kineticfact = size/hubble/rng/amax*a*100.E5*hubble*pow(amax/a,1.5); */
	/* 10.24*h^-1/(amax/a)*a*100km/s/Mpc*h*(1+z)^1.5 
	r2kineticfact = size/hubble/amax*a*a*100.E5*hubble*pow(amax/a,1.5); */
	/*
	 * X : in real scale, dX/dt : peculiar velocity in real scale
	 * dx : displacement in simulation dimension.
	 * da : simulation time step
	 * dx/da : velocity in simulation dimensions
	 * dX/dt = \dot{a} dX/da
	 *       \dot{a} = a \times H_o \times sqrt{Omega(1+z)^3+...}
	 *       dX/da = Lsize*a/amax dx/da
	 */

	com2real = 1.L/h0 * a/amax*Mpc;
	real2com = 1.L/com2real;
//	printf("pntmass = %g\n",pntmass);
	pntmass = Msun/h0;
	potentfact = G*pntmass/com2real;
//	r1kineticfact = a/amax*100.E5L*Hsub;
	r1kineticfact = 100.E5L*Hsub *a/amax;
	r2kineticfact = 1.E5L;


	if(myid ==0){
		fprintf(stdout,"P%d has omep=%g rng=%g r1kineticfact=%g r2kineticfact=%g potentfact=%g\n",myid,omep,rng,r1kineticfact,r2kineticfact,potentfact);
		fflush(stdout);
	}
	return 0;
}
