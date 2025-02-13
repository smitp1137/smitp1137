
#include<stdio.h>
#include<math.h>
#include<dir.h>
#include<sstream>
#include<iostream>
#include<fstream>

using namespace std;

//DEFINING CONSTANTS
#define NI 102
#define NJ 102
#define pi 3.141592
#define L 1
#define H 1
#define DELT 0.001

//DEFINING FUNCTION
void grid_generation();
void initial_condition();
void boundary_condition();
void coefficients_calc();
void solve_implicit();
void WRITE_FILE();
void update();

//DEFINING GLOBAL VARIABLES

int i, j;
double alpha=0.01, rho=1, cp=1,q_gen=0;
double X0[NJ][NI], X1[NJ][NI], X2[NJ][NI], X3[NJ][NI], Y0[NJ][NI], Y1[NJ][NI], Y2[NJ][NI], Y3[NJ][NI], XCELL_X[NJ][NI], XCELL_Y[NJ][NI]; // Nodes and Vertices
double s0[NJ][NI], s1[NJ][NI], s2[NJ][NI], s3[NJ][NI], DXF0[NJ][NI], DXF1[NJ][NI], DXF2[NJ][NI], DXF3[NJ][NI], VOL[NJ][NI]; // Surface Area and Nodal Distance
double DELX = L / (NI-2.0),DELY = H / (NJ-2.0);
int NCELLI = NI-1, NCELLJ = NJ-1,NSTEP=0;
double PHI[NJ][NI], PHI_OLD[NJ][NI], PHI_INT[NJ][NI]; // Temperature Array
double ae[NJ][NI], as[NJ][NI], aw[NJ][NI], an[NJ][NI], ap[NJ][NI]; // Coefficients

int main()
{
	mkdir("SMIT");
	double TMAX=5.01;
	int WRITE_INT_CONTOUR=10;
	grid_generation();
	initial_condition();
	boundary_condition();
	WRITE_FILE();
	
	coefficients_calc();
	for(NSTEP=1; NSTEP<=(TMAX/DELT); NSTEP++)
	{
		boundary_condition();
		solve_implicit();
		if( NSTEP % WRITE_INT_CONTOUR == 0)
		{
			WRITE_FILE();
		}
		
		update();
	} 
	return 0;
}

	void grid_generation() //defining the coordinates of interior nodes with central nodes
	{
		for (j = 1;j < NCELLJ;j++)
	{
	
		for (i=1;i < NCELLI;i++)
		{
		    X0[j][i] = DELX * (i - 1);
			X1[j][i] = DELX * i;
			X2[j][i] = DELX * i;
			X3[j][i] = DELX *(i-1);
			
			Y0[j][i] = DELY * (j - 1);
			Y1[j][i] = DELY * j;
			Y2[j][i] = DELY * j;
			Y3[j][i] = DELY * (j-1);

			XCELL_X[j][i] = 0.25 * (X0[j][i] + X1[j][i] + X2[j][i] + X3[j][i]);
			XCELL_Y[j][i] = 0.25 * (Y0[j][i] + Y1[j][i] + Y2[j][i] + Y3[j][i]);
		}
	}
	 
	//defining the coordinates of boundary nodes: Counter Clockwise (W S E N)
	
	for (j=1;j < NCELLJ;j++)
	{
		//WEST WALL 
		XCELL_X[j][0] = 0;
		XCELL_Y[j][0] = XCELL_Y[j][1];
		
		//EAST WALL
		XCELL_X[j][NCELLI]=L;
		XCELL_Y[j][NCELLI] = XCELL_Y[j][NCELLI - 1];
	}
	
	for (i = 1;i < NCELLI;i++)
	{
		//NORTH WALL
		XCELL_X[NCELLJ][i] = XCELL_X[NCELLJ - 1][i];
		XCELL_Y[NCELLJ][i] = H;
		
		//SOUTH WALL
		XCELL_X[0][i] = XCELL_X[1][i];
		XCELL_Y[0][i] = 0;
	}
	
	//defining the surface area and nodal distance
	
	for (j = 1;j < NCELLJ;j++)
	{
		for (i = 1;i < NCELLI;i++)
		{
			s0[j][i] = -DELY;
			s1[j][i] = -DELX;
			s2[j][i] = DELY;
			s3[j][i] = DELX;

			DXF0[j][i] = DELX;
			DXF1[j][i] = DELY;
			DXF2[j][i] = DELX;
			DXF3[j][i] = DELY;

			VOL[j][i] = DELX * DELY;
		}
	}
}
	
	void initial_condition()
	{
	
		for (j = 0;j <= NCELLJ;j++)
		{
			for (i = 0;i <= NCELLI;i++)
			{
				PHI_OLD[j][i] = 100;
				PHI[j][i] = 100;
			}
		}
	}
	
	void boundary_condition() //Considering Dirichlet Boundary Condition on every wall
	{
	for (j = 1;j < NCELLJ;j++)
	{
		// West wall
		PHI[j][0] = 300; 
		PHI_OLD[j][0] = 300;

		// East wall
		PHI[j][NCELLI] = 100; 
		PHI_OLD[j][NCELLI] = 100;
	}
	
	for (i = 1;i < NCELLI;i++)
	{
		// South wall
		PHI[0][i] = 200; 
		PHI_OLD[0][i] = 200;

		// North wall
		PHI[NCELLJ][i] = 150;
		PHI_OLD[NCELLJ][i] = 150;
	}
} 

        void coefficients_calc()
        {
			for (j = 1;j < NCELLJ;j++)
	        {
		for (i = 1;i < NCELLI;i++)
		{
			aw[j][i] = fabs((alpha * s0[j][i]) / DXF0[j][i]);
			as[j][i] = fabs((alpha * s1[j][i]) / DXF1[j][i]);
			ae[j][i] = fabs((alpha * s2[j][i]) / DXF2[j][i]);
			an[j][i] = fabs((alpha * s3[j][i]) / DXF3[j][i]);

			ap[j][i] = aw[j][i] + as[j][i] + ae[j][i] + an[j][i];

		}
	}
}

		void solve_implicit ()
		{  
           

                    for (j = 0;j <= NCELLJ;j++)
                {
                    for (i = 0;i <= NCELLI;i++)
                    {
                        PHI_INT[j][i]=PHI[j][i];

                    }
                }
                  
					float l=1;
					float error=0;
                   while (l>0.01)
                   {		l=0;
                            for (j = 1;j < NCELLJ;j++)
                        {
                            for(i=1;i < NCELLI;i++)
                            {
                                PHI[j][i] = ((q_gen*DELT)/(rho*cp)) + PHI_OLD[j][i] + (DELT/VOL[j][i])*(aw[j][i]*PHI_INT[j][i-1] + ae[j][i]*PHI_INT[j][i+1] + as[j][i]*PHI_INT[j-1][i] + an[j][i]*PHI_INT[j+1][i] - ap[j][i]*PHI_INT[j][i]);
                                l=error+(pow((PHI[j][i]-PHI_INT[j][i]),2));
                                l=sqrt(error/10000);
								PHI_INT[j][i]=PHI[j][i];
                            }
                        }
                   }
                   

                   update();
               
            
            }
			
			
		
		
		void update()
		{
	for (j = 1;j < NCELLJ;j++)
	{
		for (i = 1;i < NCELLI;i++)
		{
			PHI_OLD[j][i] = PHI[j][i];
		}
	}
}

		void WRITE_FILE()
		{
			int i,j;
			ofstream file;
			stringstream ss;
			ss << "SMIT/TEMP_PROFILE_TIME_" << NSTEP*DELT << ".dat";
			file.open(ss.str().c_str());
			file<<"VARIABLES = "<<'"'<<"X"<<'"'<<", "<<'"'<<"Y"<<'"'<<", "<<'"'<<"TEMPERATURE"<<'"'<<endl;
			file<<"ZONE I="<<NI<<", J="<<NJ<<", F=POINT"<< endl;
			file<<"SOLUTIONTIME="<<NSTEP*DELT<< endl;
			
			for(j=0;j<=NCELLJ;j++)
			{
				for(i=0;i<=NCELLI;i++)
				{
			file<<XCELL_X[j][i]<<" "<<XCELL_Y[j][i]<<" "<<PHI[j][i]<<endl;
		}
	}
}
    
	
 



