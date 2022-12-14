//********************* 2D Diffusion Solver *****************
//Vishal Indivar Kandala | 6/25/21
//**********************************************************
//Finite Difference method: Central Differencing is implemented
//***********************************************************
// DMDA Data structure available in PETSc is utilized.
// *********************************************************
//Domain: L =0.4; W=0.3	    |Stencil: Star
//		____T3____  |
//	T1=313	|    |	 |  |	 	 T(i,j+1)
//	T2=273	|____|_W_|  |	    	    |			
//	T3=283	T2   L	 T4 |	    	    |
//	T4=273	|    |	 |  | T(i-1,j)--- T(i,j)---T(i+1,j)
//		|___T1___|  |		    |
//			    |		    |
//			    |		 T(i,j-1)    
// ******************************************************************
// 			 T_xx+T_yy=0
// *****************************************************************
// 		2*(F_x+F_y)*T(i,j)-(F_x*T(i+1,j))-(F_x*T(i-1,j))-(F_y*T(j,j+1))-(F_x*T(i,j-1))
// *****************************************************************
// *****************************************************************
//
#include<petsc.h>
#include<stdlib.h>
#include<stdio.h>
//************************************************************************
//Function to create the RHS vector b
PetscErrorCode FormRHS(DM da,Vec b)
{
	int i,j;
	
	double hx, hy, BC[4]={313.0,273.0,283.0,273.0}, one = 1.0;

	double **ab;
	
	DMDALocalInfo info;

	DMDAGetLocalInfo(da,&info);

	DMDAVecGetArray(da,b,&ab);

	hx=0.4/(info.mx-1);
	
	hy=0.3/(info.my-1);

	

	//BC[0]=313.0; BC[1]=273.0; BC[2]=283.0; BC[3]=273.0;
	
	for  (j=info.ys;j<(info.ys+info.ym);j++)
	{	
		
		for (i=info.xs;i<(info.xs+info.xm);i++)

		{
			
			ab[j][i]=0;
			
			if(i == 0) 
			{	
				if(j == 0)
				{
					ab[j][i]=(BC[0]+BC[1])/2; // Numerical Decoration.
				}
				if(j == info.my-1)
				{
					ab[j][i]=(BC[0]+BC[2])/2;
				}
				if (j>0 && j<info.my-1)
				{
					ab[j][i]=BC[1];
				}
			}
			if (i == 1)
			{
				ab[j][i]=ab[j][i]+BC[1]/(hx*hx);
				if(j == 0)
				{
					ab[j][i]=BC[0];
				}
				if ( j == 1)
				{
					ab[j][i]=ab[j][i]+(BC[0]/(hy*hy));
				}
				if (j == info.my-2)
				{
					ab[j][i]=ab[j][i]+(BC[2]/(hy*hy));
				}
				if (j == info.my-1)
				{
					ab[j][i]=BC[2];
				}
			}
			if (i>1 && i<info.mx-2)
			{
				if(j == 0)
				{
					ab[j][i]=BC[0];
				}
				if (j == info.my-1)
				{
					ab[j][i]=BC[2];
				}
				if ( j == 1)
				{
					ab[j][i]=ab[j][i]+BC[0]/(hy*hy);
				}	
				if (j == info.my-2)
				{
					ab[j][i]=ab[j][i]+BC[2]/(hy*hy);
				}
				if (j>1 && j<info.my-2)
				{
					ab[j][i]=0;
				}
			}
			if (i == info.mx-2)
			{
				ab[j][i]=BC[3]/(hx*hx);

				if(j == 0)
				{
					ab[j][i]=BC[0];
				}
				if ( j == info.my-1)
				{
					ab[j][i]=BC[2];
				}
				if ( j == 1)
				{
					ab[j][i]=ab[j][i]+(BC[0]/(hy*hy));
				}	
				if (j == info.my-2)
				{
					ab[j][i]=ab[j][i]+(BC[2]/(hy*hy));
				}
			}	
			if (i == info.mx-1)
			{
				if(j == 0)
				{
					ab[j][i]=ab[j][i]+(BC[3]+BC[0])/2;
				}
				if(j == info.my-1)
				{
				
					ab[j][i]=ab[j][i]+(BC[3]+BC[2])/2;
				}
				if (j> 0 && j< info.my-1)
				{
					ab[j][i]=ab[j][i]+BC[3];
				}
			}	 // End of If-Else clause	
	
		} //End of y-loop
	}//End of x-loop
	
	DMDAVecRestoreArray(da,b,&ab);
	
	return 0;
}
//***********************************************************************
//Function to Create the matrix A
//************************************************************************
PetscErrorCode FormMat(DM da, Mat A)
{
	DMDALocalInfo info;
	
	int i,j,ncols;
	
	double hx,hy,v[5];

	MatStencil row,col[5]; // links the da and matrix A

	DMDAGetLocalInfo(da,&info);

	hx=0.4/(info.mx-1);
	
	hy=0.3/(info.my-1);

	// The (m*n)*(m*n) A matrix is solved to find the rasterized (m*n)*(1) vector solution in PETSc as the solution cannot be stored as a matrix in PETSc.

	for (i=info.xs;i<(info.xs+info.xm);i++)
	{	
		for (j=info.ys;j<(info.ys+info.ym);j++)
		{
		
			row.j = j; row.i = i; // The i*j th row in the (m*n)*(m*n) A matrix is selected.
			
// No.of entries to be made in each row of A matrix	col[0].j = j; col[0].i =i; // Using the matrix Stencil Structure, we can place the values in v at the location (row.i*row.j),(col.i*col.j) so we are selecting the indices and placing each value.
			ncols = 1;
 
			col[0].i=i; col[0].j=j;
			
			v[0] = 2*((1/(hx*hx))+(1/(hy*hy)));
				
			// Addressing Boundaries

			if(i>1)
			{
				col[ncols].i =i-1; col[ncols].j = j;
				v[ncols] = -1/(hx*hx);
				ncols=ncols+1;
			}
			if(j>1)
			{	
				col[ncols].i=i; col[ncols].j=j-1;
				v[ncols]=-1/(hy*hy);
				ncols=ncols+1;
			}
			if(i<info.mx-2)
			{
				col[ncols].i=i+1; col[ncols].j=j;
				v[ncols]=-1/(hx*hx);
				ncols=ncols+1;
			}
			if(j<info.my-2)
			{
				col[ncols].i=i;col[ncols].j=j+1;
				v[ncols]=-1/(hy*hy);
				ncols=ncols+1;
			}
			if (i == 0 || j == 0 || i == info.mx-1 || j == info.my-1)
			{
				ncols=1;
				v[0]=1;
	
			} // End of if-else clause
		

			MatSetValuesStencil(A,1,&row,ncols,col,v,INSERT_VALUES); // Enter the colums (the entire row) into the matrix A associated with da.	
			
		} // End of y-loop
		
	} // End of x-loop

	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);

	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

	return 0;

}
//************************************************************
//Function to write the solution to a CSV file
//***********************************************************
PetscErrorCode WriteSol(DM da, Vec u)
{
	int i;
	
	int j;

	double **au,v;

	FILE *fid;

	const char * dir = "../data/";
	
	const char * fname = "2diff";

	const char * ftype = ".csv";
	
	char name_buffer[4096];
	
	DMDALocalInfo info;
	
	DMDAGetLocalInfo(da, &info);

	DMDAVecGetArray(da,u,&au);
	
	sprintf(name_buffer,"%s%s-%dx%d%s",dir,fname,(int)info.mx,(int)info.my,ftype);	

	fid=fopen(name_buffer,"w+");

	for(j=info.ys;j<(info.ys+info.ym);j++)
	{
		for(i=info.xs;i<(info.xs+info.xm);i++)
		{
			v=au[j][i];
			if(i==info.mx-1)
			{
				fprintf(fid,"%.2lf",v);
			}	
			else
			{
				fprintf(fid,"%.2lf,",v);
			}
		}
		fprintf(fid,"\n");
	}
	
	DMDAVecRestoreArray(da,u,&au);

	fclose(fid);


	return 0;
}

//************************************************************
//Main
//************************************************************

int main(int argc, char **argv)
{	
	int flg1,flg2,flg3,its,n,m,p1,p2,p3;

	const char * ksp_type;
	
	PetscErrorCode ierr;

	PetscLogDouble duration,start, stop;

	DM da; //Distributed Member (or ) Data Management object. (Used for a structured grid).

	Mat A;

	FILE *fidtime;

	const char * dir = "../data/";

	const char * timefilegs = "timegs";

	const char * timefilegmres = "timegmres";

	const char * timeftype = ".csv";

	char name_buffer[1024];

	Vec b,u;

//	Vec useq;

//	VecScatter ctx;	

	KSP ksp;

	DMDALocalInfo info; // Data Structure to hold information about the array
	
	//n=atoi(argv[1]); m=atoi(argv[2]);
	
	PetscInitialize(&argc,&argv,NULL,"Solve Poisson Equation in 2D");

	ierr=PetscTime(&start); CHKERRQ(ierr);
	//***************************************************
	// Creating the DMDA Data Structure and declaring it over the matrix A
	//****************************************************
	ierr=DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,2000,2000,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&da);CHKERRQ(ierr);

	ierr=DMSetFromOptions(da);CHKERRQ(ierr);
	
	ierr=DMSetUp(da);CHKERRQ(ierr);

	ierr=DMCreateMatrix(da,&A); CHKERRQ(ierr);// The size of A is determined by grid size specified in DMDACreate
	
	ierr=MatSetFromOptions(A);CHKERRQ(ierr);

	//**************************************************
	// Creating the vectors b and u from DMDA object
	// *************************************************

	ierr=DMCreateGlobalVector(da,&b);CHKERRQ(ierr);

	ierr=VecSetFromOptions(b);CHKERRQ(ierr);

	ierr=VecDuplicate(b,&u); CHKERRQ(ierr);// Initializing u

	//*********************************************
	//Form Matrix A
	//*********************************************
	ierr=FormMat(da,A); CHKERRQ(ierr);

	//*********************************************
	//Form vector b
	//********************************************
	ierr=FormRHS(da,b); CHKERRQ(ierr);
	
//	PetscPrintf(PETSC_COMM_WORLD, "***** b *********** \n");

//	VecView(b,PETSC_VIEWER_STDOUT_WORLD);
	//*******************************************
	//Solve the system Ax=b
	//******************************************
	ierr=KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);

	ierr=KSPSetOperators(ksp,A,A); CHKERRQ(ierr);

	ierr=KSPSetFromOptions(ksp); CHKERRQ(ierr);

	ierr=DMDAGetLocalInfo(da,&info);	
	
	ierr=KSPSolve(ksp,b,u); CHKERRQ(ierr);
	
	ierr=KSPGetTotalIterations(ksp,&its); CHKERRQ(ierr);
	
	ierr=KSPGetType(ksp,&ksp_type);

	ierr=PetscTime(&stop); CHKERRQ(ierr);

	duration=stop-start;

	flg1=strcmp(ksp_type,"richardson");
	
	flg2=strcmp(ksp_type,"preonly");
	
	flg3=strcmp(ksp_type,"gmres");
	
	if(flg1==0 || flg2==0)
	{
		sprintf(name_buffer,"%s%s%s",dir,timefilegs,timeftype);
	}
	else if(flg3==0)
	{
		sprintf(name_buffer,"%s%s%s",dir,timefilegmres,timeftype);
	}
	
	fidtime=fopen(name_buffer,"a");	

	ierr=DMDAGetInfo(da,NULL,NULL,NULL,NULL,&p2,&p1,&p3,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
	
	ierr=PetscFPrintf(PETSC_COMM_WORLD,fidtime,"%d,%d,%d,%d,%.2lf\n",its,p2+p1+p3,info.mx,info.my,duration);

//	************************************************
//	Print the solution out to a file 
//	************************************************

//	ierr=VecScatterCreateToZero(u,&ctx,&useq); CHKERRQ(ierr);

//	ierr=VecScatterBegin(ctx,u,useq,INSERT_VALUES,SCATTER_FORWARD);

//	ierr=VecScatterEnd(ctx,u,useq,INSERT_VALUES,SCATTER_FORWARD);
	
	ierr=WriteSol(da,u); CHKERRQ(ierr);
	
	//******************************************
	// Memory De-allocation for vectors u,uexact,b, matrix A, solver object ksp and DM Data Structure da
	// ****************************************
	ierr=VecDestroy(&u);CHKERRQ(ierr);

//	ierr=VecScatterDestroy(&ctx);CHKERRQ(ierr);

//	ierr=VecDestroy(&useq);CHKERRQ(ierr);

//	ierr=VecDestroy(&useq);CHKERRQ(ierr);

	ierr=VecDestroy(&b);CHKERRQ(ierr);

	ierr=MatDestroy(&A);CHKERRQ(ierr);

	ierr=KSPDestroy(&ksp);CHKERRQ(ierr);

	ierr=DMDestroy(&da);CHKERRQ(ierr);

	return PetscFinalize();
}
