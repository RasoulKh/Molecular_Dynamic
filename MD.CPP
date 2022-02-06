// Added a coefficient to the velocity to apply the thermostate
// Same as MS_MD18.cpp but changed location of "New_Cells" To improve Velocity"  (19)
// Redefinition of Temperature (19_1)
// addition of a velocity dependent term to the force	(19-2)
//	change boundary check in integration functions (B&V) (19-2)
// improve speed of peridic boundary by using "else"     (19-3)
// thermal wall (19-3)
//	corrected1 in MAKE_Neighb_Cell()   (19-3)
//  Computation of properties in cells	(19-4) by adding : Sampling(),Sam_Cell_Det(int i1),Make_Sam_Cell_Member(),Devide_Geo() (19-4)
//	delete coefficient to the velocity to apply the thermostate (19-5)
//	add Write_Sampled_Data() (19-5)
//	add Exreral_Force()    (19-5)
//  two types of neighbor list generation : for small systems, All_Pair and for large systems using cell subdivision		(19-5)
//	Avoid reconsideration of considered cells in making neighbor lists	(19-6)
//  N_Sam_CELL[0] >> N_Sam_CELL[1] >> N_Sam_CELL[2] in Input()           (19-6)
//  added density number to sampling()    (19-8)
//	adding an if to Sampling ()  so it can hande finer mesh   (19-9)
//	Extended Compute_Force() to Van der Waals, Bond, Angle, Dihedral, Impropers, and Electrostatic forces   (19-10)
//  Read Sampling_start  & Sampling_step from inout file   (19-10)
//	correction of temperature redef. by deviding v_mean by n_mesure and commenting line 674(19-11)
//	deleted additional "for" in calc_prop()      (19-11) 
//	some terms (Pot_tot=0,..., Force=0 & external_force() moved fro LJ_force() to Compute_force()	(19-11)
//	in the Single_Step() add the condition (Thermp_type=S) before scaling command     (19-12)
//	move Show_prop() and calc_properties from integration functions to the single_step()	(19-12)
// remove a bug in Kin_tot_Measure+=...         (19-12)
// add   Kin_tot+=  to Calc_Properties()      (19-12)
//	add void make_Inf_Solid_list() & void Inf_Solid_Body()  to functions and Inf_Solid_Body() to Calc_Force()	(19-13)
//  corrected V_mean calculation  (19-4)
//  definition Immobalize_system() & adding it to Initializing()    (19-14)
//  Adding Epsilon_MAT[][] & Sigma_MAT[][] which reduces computational cost     (19-15)
//	applying 3rd Newton law concept in neighbor list construction          (19-15-1)
//  using 2*delr_max concept for reconstruction of neighbor lists           (19-15-2)
//	adding a new loop for Inc determination in LJ_Force() due to Vel_Set            (19-15-2)
//	adding a new loop for Inc determination in LJ_Force()  due to Inf_Solid_Body           (19-15-2)
//	aeeing Write_Trajectory_Restart()and read_Restart();                                          (19-15-3)
// Fixing a bug in void	Create_Neighbourhood();               
//					if (k4!=Cell[i]	&& (Considered_Cell[k4]==0)) {
//				Considered_Cell[k4]=1;       (19-15-3)
//	Separation of data file and iput file  (19-15-5)
//	Including Bond and Angle		(19-15-5)
//	Fixing a bub in Input of angle type parameters  (19-15-5)
//  adding bond and angle Energies to Potential Energy   (19-15-5)
//	some changes in LJ_"Force() or excluding atoms  (19-15-7)
//  Introducing Diplay_Scale                    (19-15-7)
//	Computation of Lacal temperature in sampling and output() (19-15-8)
//	redefinition of local tamp     (19-15-9)
//	delete the condtition (if Cell[i]==Cell[j] else) from neighbor list construction.  (19-19-10)
//	excluding some of cells which those members are surely in the neighbor list.  (19-15-10)  **
//  change search from atom oriention to cell oriention   (19-15-11)
//	remove **   (19-15-11)
//  excluding od pairs in the same cell, angle, inf body and constant vel is moved from LJ_force() to create_Neighberhood()  (19-15-12-1)
//  checking for rij < r_cut is added to LJ  (19-15-12-1) 
//  added Histogram    (19-15-12-2)
//						(19-15-12-2-1)
//  deleting "while" from peridic boundary condition (19-15-12-2-2)
// adde Xi to the output file (19-15-12-2-2)
// commented check for scsling velocity in single_step  (19-15-12-2-2)
// commented unnecessary fix_velocity in Integrate_V  (19-15-12-2-2)
// separate in and out H2 in histogram   (19-15-12-2-5)
// added RDF  (19-15-12-2-7)
// decrease memory usage    (19-15-12-2-7)
// some improvments in LJ    (19-15-12-2-7)
// in make_neighbor2 the last sentence is commented  (19-15-12-2-7)
// in the last line of Allocate_arrays which we used nim_of_atoms/4, in can make a bug for systems with small side walls. 
// some changes in sampling()    (19-15-12-2-8)
// some changes in sampling() in temperature , addition N_T_S    (19-15-12-2-9)
// some changes in sampling()  , deleting N_T_S    (19-15-12-2-10)
// replace int N_atoms_passed... by double M_atoms_passed...  (19-15-12-2-10)
// Unified Codes  (20)
// Rebuilt duration is deleted (21)
// n is deleted for the inf_solid in input file  (21)
//////////////////////////////////////////////////////////
//	3D
//  24_4:  a single atom Neighbor list making is developed behind the previous overal Neighbor list making
//  24_5:  Refresh_Cell() is improved
//	24_5:  Improve P.B.C in Integrate_V
//  24-6   Dynamic determination of number of atoms in the cells and sam cells
//	24-6.1 the code improved by Dihedral inter_atomic forces and Thrus force
//  24-6-2 Inulated Wall added to the Code and related changes     5/25/2010
//  24-6-3 Inf_solid -----> Solid and two related functions changed  may 30 2010
//  24-6-3 calculation of Max_N_Of_Sam_Cell changed    may 31 2010
//  24_6_3 max_N_of_cell calculated and added as 24_6   june 16 2010
//  24_6_4 Variable Atom Number  July  4  2010
//  24_6_5 GetInputs -------> Get_Input+Get_Data  July  4  2010
//  24_6_5 Some aalocation moved from allocatearrays to Get_Input    July  4  2010
//  24_6_5_01 Num_Of_Atoms added to Read and Write_restarts    July  5  2010
//  24_6_5_02	Deleting Bonds and Angles     added to DeleteAtoms         July  5  2010
//  24_6_5_02	X_Low,Y_Low,Z_Low -------->     X_Low_E,Y_Low_E,Z_Low_E in  sam_cell_det         July  5  2010
//  24_6_5_02	atomtype  added to read and write restart       July  5  2010
//	24_6_5_02 X_Low,Y_Low,Z_Low , High -------->     X_Low_E,Y_Low_E,Z_Low_E ad High_E in  sam_Cell_Det  and Devide_Geo       July  5  2010
//  24_6_5_03  Rho_max reduced to 8  July  6  2010
//  24_6_5_03  sigma_min and sigma_max determination moved from make_naighb_cell to get_input  July  6  2010
//   24_6_5_03          max_ng  is globalized       July  6  2010
//   24_6_5_03    Refresh_Del_list added to integrate_V    July  6  2010
//   24_6_5_03    an if removed from Dump_trajectory    July  6  2010
//   24_6_5_03    an = -----> == in an if(  )    July  6  2010
//	 24_6_5_6 Inlets added to the code 	July  10  2010
//	24_6_5_6   X_Low ----> X_Low_E in  X_Of_Sam_Cell  and so for Y_Of_Sam_Cell and  Z_Of_Sam_Cell      July  12  2010
//     deleting and NULLing for Inlet_Atoms  in Inlets()          July  14  2010
//   24_6_5_6_1   pot devided between i and j in LJ_Force()           July  14  2010
//	 24_6_5_6_4 	Normal_Random(double, double) added to code             July  19  2010
//     24_6_5_6_4          some  corrections for Pressure_work in LJ and LJ_of_Atom               July  20  2010
//     24_6_5_6_4          last line of   LJ_of_Atom is multiplied by 2             July  20  2010
//    24_6_5_6_5       X_Low_E=X_Low   and so on.                  July  21  2010
//    24_6_5_6_5       if (Atoms[m].tp==tp) added for temp and vel setting in Inlets()                  July  21  2010
//    24_6_5_6_5          Inlet_Atoms deleted    in Inlets()                  July  21  2010
//    24_6_5_7          Outlet added                 July  21  2010

//    24_6_5_7_2                 v_mean[0]=0;   and so removed from calc_properties             Agust 10 2010
//    24_6_5_7_2         Calc_Properties_Mid() and	Calc_Properties() added before immoblize        Agust 10 2010
//       24_6_5_7_2                Velocity_Profile()   deleted          Agust 10 2010
//      24_6_5_7_2      Create_Neighbourhood2     deleted              Agust 10 2010
//     24_6_5_7_2 		numofNeighbors localized           Agust 16 2010
//        24_6_5_7_2           for (i2=1;i2<Max_N_Of_Sam_Cell;i2++) Sam_CELL_MEMBER[i1][i2]=-1;    commented                Agust 19 2010
//           24_6_5_7_2                 deinition of Max_N_Of_Cell  changed                 Agust 19 2010
//            24_6_5_7_2             Refresh_Delte List Moved               Agust 25 2010
//              24_6_5_7_2       Harmonic_Angle_Parameter[i].K*=(180*180)/(pi_N*pi_N);      aded            Agust 30 2010
//              24_6_5_7_2      Some "if"s changed to "else if" for priodic B.C.           Agust 30 2010
 //              24_6_5_7_2    some therostat parameteres deleted from input file      Agust 30 2010
//       24_6_5_7_2        Trajectory_Write_Step Properties_Write_Step Restart_Write_Step added        Sept 2 2010
//         24_6_5_7_2        for (i2=0;i2<Outlet[i1].N_C;i2++) ---->   for (i2=0;i2<Inlet[i1].N_C;i2++)        Sept 4 2010
//       24_6_5_7_2    Max_N_Of_Sam_Cell  calculation corrected        Sept 6 2010
//////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <string.h>
#include <direct.h>

using namespace std;
//////////////////////////////////////
//CONSTANTS //
double	Boltzmann_const=8.314e-7,coulomb_factor = 1.0, dielectric_constant = 138.7008e-3;
int MMM=3;
double Rho_max=8.0;
double pi_N=4*atan(1);
int max_ng=4100; 
// Structures 
//////////////////////////////////////
struct C_T_W {
    int tp;   
	double vx;
	double vy;
	double vz;
	double Temp;
	double K;
	double xi;
 };
//////////////////////////////////////
struct C_F_W {        
    int tp;   
	double vx;
	double vy;
	double vz;
	double Flux;
	double K;
 };
/////////////////////////////////////
struct F_Struct {
    int tp;   
    char dr;   
    double  mag;     
 };
/////////////////////////////////////
struct Bound_T {
	char x;
	char y;
	char z;
};
////////////////////////////////////
struct V_Set {
	int tp;
	double vx;
	double vy;
	double vz;
};
////////////////////////////////////
struct MORSE_B {
	double D;
	double r0;
	double alpha;
};
///////////////////////////////////
struct Harmonic_B {
	double K;
	double r0;
};
///////////////////////////////////
struct DIHEDRAL {
	double K;
	int multiplicity;
	double cos_shift;
	double sin_shift;
	int sign;
};
///////////////////////////////////
struct C_S_A {
	double K;
	double theta0;
};
///////////////////////////////////
struct Inlet_Struct {
	int atom_tp;
	double Rho;
	double temp;
	double Vx;
	double Vy;
	double Vz;
	double x1;
	double y1;
	double z1;
	double z2;
	double x2;
	double y2;
	int N_C;    // Number of cells belonging to 
	int Cells[1000];
	double Cell_N_Res[1000];
};
///////////////////////////////////
struct Outlet_Struct {
	char dir;
	double x1;
	double y1;
	double z1;
	double z2;
	double x2;
	double y2;
	int N_C;    // Number of cells belonging to 
	int Cells[100];
};
///////////////////////////////////
struct Solid {
	int Num;
	int tp;
	double init_Vx;
	double init_Vy;
	double init_Vz;
	double fx;
	double fy;
	double fz;
	double si;
	double I;
};
///////////////////////////////////
struct Propeller {
	int Atom_No;
	int Atom_No_Direction_1;
	int Atom_No_Direction_2;
	double fx;
	double fy;
	double fz;
	double F;
	double xx;
	double yy;
	double zz;
 };
///////////////////////////////////
struct Atom_Struct {
	double	R_Vec[3];
	double	V_Vec[3];
	double	A_Vec[3];
	double	force[3];
	double	Pot;
	double  Pres_Work;	
	int		MoleculeOfAtom;
	int		Num_Of_Bonds_Of_Atom;
	int		Num_Of_Angles_Of_Atom;
	int		Num_Of_Dihedrals_Of_Atom;
	int		*Bonds_Of_Atom;
	int		*Angles_Of_Atom;
	int		*Dihedrals_Of_Atom;
	int		tp;
	int		*NeighborList;
	int		Cell;
	int		Sam_Cell;
	int		Delete_List;
	int	    Inc_vec;    //used for inclusion of atoms in calc_prop
};
/***********************/
/***	VARIABLES    ***/
/***********************/
////////////////////
// FOR SIMULATION //
////////////////////
double	Time,DeltaT;
int		MaxStep,Trj_Write_Step, Prop_Write_Step, Restart_Write_Step,StepNum=0;
int		Num_Of_Atoms,Num_Of_Atom_Types,Num_Of_Bonds,Num_Of_Bond_Types,Num_Of_Angles,Num_Of_Angle_Types,Num_Of_Dihedrals,Num_Of_Dihedral_Types;
int		Num_Of_Thrust_Sets,Num_Of_Force_Sets,Num_Of_Solid_Sets;
int		N_freezed_mols,Num_Of_Vel_Sets,*freezed_mols=NULL;
int		**Solid_Atoms=NULL;
bool	ContinueCondition,ReCalc_Neighb=1; 
double	R_Cut_Off,Delta_R,R_Cut_Off_R,Delta_R_R ;
double  X_Low,Y_Low,Z_Low,X_High,Y_High,Z_High,X_H_L,Y_H_L,Z_H_L;
double	Pot_tot,Kin_tot,v_mean[3],V_mean=0.; 
double  M_system=0.0;
double  W_SR,Sys_Temperature,Sys_Pressure,Volume,Target_Temp,Xi=0.0,Xip=0.0,T_Coupling;
char	Thermo_tp,Wall_Thermo_Method='S';
char	Integ_Method='V';
double	rebuild_neighb_char=0,sigma_min,sigma_max;
char	Trj_INIT_MODE,Smp_INIT_MODE;
double  Display_Scale;
double  Dx[3];
double	X_L_Measure,X_H_Measure,Y_L_Measure,Y_H_Measure,Z_L_Measure,Z_H_Measure;
int		N_Measure;
int		Dimension=3,N_CELL[3],N_CELLS,N_Sam_CELLS,N_Sam_CELL[3],Sampling_Step,Sampling_start;
double  Dx_s[3];
int	    Sam_Step=0;
double  R_Cut_Coul;
double  max_R_Cut, Rho_Water;
int     Water_Mol=0,Coulomb_OK;
//    for Variable atom number
int		Max_Num_Of_Atoms;
double  Mem_Increase_factor;
int		Delete_Bound_OK;
double  X_Low_E,Y_Low_E,Z_Low_E,X_High_E,Y_High_E,Z_High_E;
///////////////////
// FOR EACH ATOM //
///////////////////
Atom_Struct	*Atoms=NULL;
////////////////////////
// FOR EACH ATOM TYPE //
////////////////////////
double	*Epsilon=NULL,*Sigma=NULL,*AtomMass=NULL,*q = NULL;
double  **Sigma_MAT=NULL,**Epsilon_MAT=NULL;
////////////////////////
//    FOR EACH Cell   //
////////////////////////
int		**CELL_MEMBER=NULL,**Neighb_Cell=NULL;
////////////////////////////
//    FOR EACH SAM CELL   //
////////////////////////////
int		**Sam_CELL_MEMBER=NULL;
double  (*Vel_Of_Sam_Cell)[3]=NULL;
double  *KE_Of_Sam_Cell=NULL,*Temp_Of_Sam_Cell=NULL;
double  *X_Of_Sam_Cell=NULL,*Y_Of_Sam_Cell=NULL,*Z_Of_Sam_Cell=NULL;
double  *M_of_Atoms_Passed_Sam_Cell=NULL;
int		*N_of_Atoms_Passed_Sam_Cell=NULL;
double	*Density_Of_Sam_Cell=NULL;
double  Volume_Of_Sam_Cell;
int		*Considered_Cell=NULL;
double  *Flux_KE_Of_X_Of_Sam_Cell=NULL,*Flux_KE_Of_Y_Of_Sam_Cell=NULL,*Flux_KE_Of_Z_Of_Sam_Cell=NULL,*Flux_pot_of_Cell_x=NULL;
double  *Flux_pot_of_Cell_y=NULL,*Flux_pot_of_Cell_z=NULL,*Flux_pres_work_Cell_x=NULL,*Flux_pres_work_Cell_y=NULL,*Flux_pres_work_Cell_z=NULL;

//////////////////////////////////////////////
//     FOR EACH BOND ANGLE AND DIHEDRAL     //
//////////////////////////////////////////////
int		(*bondlist)[3]=NULL,(*anglelist)[4]=NULL,(*dihedrallist)[5]=NULL;

//////////////////////////////////////////////
int		Num_Of_Cons_Temp_Sets,max_Wall_mol,**Cons_Temp_Atoms=NULL;
int		Num_Of_Const_Flux_Sets,**Const_Flux_Atoms=NULL;  
double  **X_e_T_Wall=NULL,**Y_e_T_Wall=NULL,**Z_e_T_Wall=NULL,*Wall_Temp=NULL,*Kin_tot_Measure_Wall=NULL;
double  **X_e_F_Wall=NULL,**Y_e_F_Wall=NULL,**Z_e_F_Wall=NULL,*Wall_Flux=NULL;  
double  *xip=NULL,*Q_heat_flux_z=NULL,*Q_heat_flux_y=NULL,*Q_heat_flux_x=NULL;
////////////////////////////////////////////////////////////////////////////////////
V_Set			*Vel_Set=NULL;
Bound_T			Boundary_Type;
F_Struct		*F_Set=NULL;
MORSE_B			*Morse_Bond_Parameter=NULL;
Harmonic_B		*Harmonic_Bond_Parameter=NULL;
Harmonic_B		*Harmonic_Angle_Parameter=NULL;
DIHEDRAL		*Harmonic_Dihedral_Parameter=NULL;
C_S_A			*Cosine_Sq_Parameter=NULL;
C_T_W           *cons_Temp_Wall=NULL;
C_F_W           *cons_Flux_Wall=NULL;          
Solid       *Solid_Wall=NULL;
Propeller		*Thrust_Set = NULL;

/////////////////////////////////////////////////////////////
int		Max_N_Of_Cell, Max_N_Of_Sam_Cell;
/////////////////////////////////////////////////////////////
// inlets
int Num_of_Inlets;
Inlet_Struct *Inlet=NULL;
/////////////////////////////////////////////////////////////
// Outlets
int Num_Of_Outlets;
Outlet_Struct *Outlet=NULL;
//*****************************************************
//*                     FUNCTIONS                     *
//*****************************************************

void	Add_Atom(double ,double,double,double,double,double,int);
void    Add_Force(int Atom_t,char direction,double Ext_Force) ;
void	Allocate_Arrays();
void	Angle_Force();
int 	Belong_Cell_to_Inlet(int,int); 
int 	Belong_Cell_to_Outlet(int,int);
void	Bond_Force();
void	Calc_Properties();
void	Calc_Properties_Mid() ;
void	Cell_Det(int i1) ;
void	Compute_Force();
void    Const_Temp_and_Flux_Wall_Force();
void	Create_Neighbourhood() ;
void	Create_Neighbourhood_of_Atom(int i);
void	Delete_Atoms() ;
void	Devide_Geo() ;
void	Dihedral_Force();
void	Dump_Trajectory();
void	External_Force();
void	Fix_Velocity() ;
void	Get_Data();
void	Get_Input();
void	Harmonic_Angle();
void	Harmonic_Bond();
void 	Harmonic_Dihedral();
void	Immobalize_System();
void	Initialization();
void 	Inlets();
void	Insert_Atom_in_Cell(int,int,double,double,double ,double ,double,double,double,double,double,double,int);
void	Integrate_V();
void	LJ_COUL_Force();
void	LJ_COUL_of_Atom(int) ;
void	Make_Cell_Member();
void    Make_Const_Temp_and_Flux_Wall_List();
void	Make_Const_Vel_List ();
void	MAKE_Neighb_Cell() ;
void	Make_Solid_List (); 
void	Make_Sam_Cell_Member();
void	New_Cells () ;
double	Normal_Random(double, double);
void	Propeller_Force();
void	Refresh_Cell(int i1,int i2,int i3);
void	Read_Sam_Restart();
void	Read__Trajectory_Restart();
void	Ref_Del_List();
void	Sam_Cell_Det(int) ;
void	Sampling() ;
void 	Set_Parameters();
void	Show_Step_Prop();
void	Single_Step();
double  Uniform_Random(double x1,double x2);
void	Usher_in_Cell(int,double,double,double ,double ,double,double ,double,double,double,double,int);
void	Van_der_Waals_Force();
void    Wall_Temp_Scaling();
void	Write_Sampled_Data();
void 	Write_Sam_Restart();
void	Write_Trajectory_Restart();
void	X_Free_Boundary(int) ;
void	X_Periodic_Boundary(int) ;
void	Y_Free_Boundary(int) ;
void	Y_Periodic_Boundary(int) ;
void	Z_Free_Boundary(int) ;
void	Z_Periodic_Boundary(int) ;
//*****************************************************
//*               INPUT & OUTPUT FILES                *
//*****************************************************
ifstream InputFile("Input.txt"); 
ifstream InputFile2("Data.txt"); 
ofstream pot_file("Potential.plt");
ofstream Drag("Drag.plt");
ofstream Tm("Time.txt");
//*****************************************************
//*****************************************************
//*                         MAIN                      *
//*****************************************************
//*****************************************************
int main ()
{  
	int starttime, stoptime, timeused;
	int ii=0;
    clock_t start=clock();
	starttime = clock();
	/***************** Potential File **********************************/
	pot_file  << "VARIABLES=" << '"' << "N" <<'"' << "," << '"' ;
	pot_file  <<"Time" << '"' << "," << '"' << "Potential_E" << '"' ;
	pot_file  <<"," << '"' << "Kinetic_E" << '"' <<"," << '"' ;
	pot_file  << "Total_E" << '"' <<"," << '"' << "Temperaure" ;
	pot_file  << '"' <<"," << '"' << "Pressure" << '"' <<"," << '"';
	pot_file  << "Xi" << '"' <<"," << '"' << "Num_ATOM"<< '"' <<'\n';
	/*******************************************************************/
	Time=0.0;
	Get_Input();
	Get_Data();
	Set_Parameters();
	Initialization();
	Make_Const_Vel_List ();
	Make_Solid_List();
	Make_Const_Temp_and_Flux_Wall_List();
	mkdir("Trajectory");
	Dump_Trajectory();
	ContinueCondition=1;
	/******************* FIRST STEP *********************/
	MAKE_Neighb_Cell();
	Devide_Geo() ;
	Make_Cell_Member();
	Create_Neighbourhood() ;
	Compute_Force();
	/****************************************************/
	while (ContinueCondition) {
		++ii;
		Single_Step();
	
		if (StepNum>=MaxStep) ContinueCondition=0;
	}
	stoptime = clock();
    timeused = stoptime - starttime;
	Tm << " T.T= " << timeused ;
	Tm.close();
	return 0;	
}
/***********************************************************************/
/***********************************************************************/
void	Single_Step()
{
	++ StepNum;
	Time+=DeltaT;
	Integrate_V();
	Inlets();
	Calc_Properties_Mid();
	Calc_Properties();
	Show_Step_Prop();
	if (StepNum%Trj_Write_Step==0) {
		Dump_Trajectory();
		if (StepNum%(Restart_Write_Step)==0) 		Write_Trajectory_Restart();
	}
	if (StepNum%Sampling_Step==0) {
		if (StepNum>Sampling_start) {
			Make_Sam_Cell_Member();
			Sampling();
		}
	}
}
/***********************************************************************/
/***********************************************************************/
void Compute_Force() 
{
	int i;
	if (rebuild_neighb_char>Delta_R_R*sigma_min/2) {
		Create_Neighbourhood() ;
		rebuild_neighb_char=0;
		cout << " Rebuilt Neighbor Lists!" << '\n'; 
	}
	Pot_tot=0;
	Kin_tot=0;
	v_mean[0]=0;
	v_mean[1]=0;
	v_mean[2]=0;
	W_SR=0;
	for (i=0;i<Num_Of_Atoms;i++) {
		Atoms[i].force[0]=0;
		Atoms[i].force[1]=0;
		Atoms[i].force[2]=0;
	}
	Van_der_Waals_Force();
	Bond_Force();
	Const_Temp_and_Flux_Wall_Force();
	Angle_Force();
	Dihedral_Force();
	External_Force();
	Propeller_Force();
}
/***********************************************************************/
/***********************************************************************/
void Van_der_Waals_Force()
{
	LJ_COUL_Force();
}
/***********************************************************************/
/***********************************************************************/
void	LJ_COUL_Force()
{
	int	i,j,jj;
	double x_ij, y_ij, z_ij, r_ij2, r_ij2_inv, r_ij6_inv, r_ij12_inv, LJ_Force=0.0, COUL_Force=0.0, Force=0.0, r_ij, x_ij2, y_ij2, z_ij2;
	double forcecoul=0.0, rinv, E_coulomb = 0.0,E_LJ=0.0, Smooth_Function=0.0, r_inv_ij2, RR_Cut_Off_R = 0.0, Charge1, Charge2;
	int  numofNeighbors;
	double sigma,epsilon;
	for (i=0;i<Num_Of_Atoms;i++) {
		Atoms[i].Pot=0; 
		Atoms[i].Pres_Work=0; 
	}
	for (i=0;i<Num_Of_Atoms;i++) {
		numofNeighbors=Atoms[i].NeighborList[0];
		for (jj=1;jj<numofNeighbors+1;jj++) {
			j=Atoms[i].NeighborList[jj];
			z_ij=(Atoms[i].R_Vec[2]-Atoms[j].R_Vec[2]);
			y_ij=(Atoms[i].R_Vec[1]-Atoms[j].R_Vec[1]);
			x_ij=(Atoms[i].R_Vec[0]-Atoms[j].R_Vec[0]);

			if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
				if ((x_ij)>X_H_L/2)  x_ij=x_ij-X_H_L;
				else if ((x_ij)<-X_H_L/2) x_ij=x_ij+X_H_L;
			}
			if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
				if ((y_ij)>Y_H_L/2) y_ij=y_ij-Y_H_L;
				else if ((y_ij)<-Y_H_L/2) y_ij=y_ij+Y_H_L;
			}
			if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
				if ((z_ij)>Z_H_L/2) z_ij=z_ij-Z_H_L;
				else if ((z_ij)<-Z_H_L/2) z_ij=z_ij+Z_H_L;
			}

			x_ij2 = x_ij*x_ij;
			y_ij2 = y_ij*y_ij;
			z_ij2 = z_ij*z_ij;

			r_ij2 = (x_ij2 + y_ij2 + z_ij2);
			r_ij = sqrt(r_ij2);
			sigma = Sigma_MAT[Atoms[i].tp-1][Atoms[j].tp-1];
			R_Cut_Off = R_Cut_Off_R*sigma;

			Charge1 = q[Atoms[i].tp-1];
			Charge2 = q[Atoms[j].tp-1];
			COUL_Force = 0.0;
			E_coulomb=0;
			if (fabs(Charge1) >0.000001) {
				if (fabs(Charge2) >0.000001) {
					if ((r_ij<R_Cut_Coul) ) {	
						r_inv_ij2 = (1.0)/r_ij2;
						rinv = sqrt(r_inv_ij2);
												forcecoul = -dielectric_constant * Charge1 * Charge2 * rinv;
						COUL_Force = (coulomb_factor*r_inv_ij2)*(-forcecoul);
						E_coulomb = coulomb_factor * (forcecoul);
					}
				}
			}

			LJ_Force = 0.0;
			E_LJ=0;
			if ((r_ij<R_Cut_Off)) {
				r_ij2_inv = (sigma*sigma)/r_ij2;
				r_ij6_inv = r_ij2_inv*r_ij2_inv*r_ij2_inv;
				r_ij12_inv = r_ij6_inv*r_ij6_inv;
				r_ij2_inv = r_ij2_inv/(sigma*sigma);
				epsilon = Epsilon_MAT[Atoms[i].tp-1][Atoms[j].tp-1];
				LJ_Force = 48*epsilon*r_ij2_inv*(r_ij12_inv-0.5*r_ij6_inv);
				E_LJ=4*epsilon*(r_ij12_inv-r_ij6_inv);
			}

 			Force = LJ_Force + COUL_Force;
		
			Atoms[i].force[0] += x_ij*Force;
			Atoms[i].force[1] += y_ij*Force;
			Atoms[i].force[2] += z_ij*Force;

			Atoms[j].force[0] -= x_ij*Force;
			Atoms[j].force[1] -= y_ij*Force;
			Atoms[j].force[2] -= z_ij*Force;
			
			
			Atoms[i].Pot+=0.5*(E_LJ+E_coulomb) ;
			Atoms[j].Pot+=0.5*(E_LJ+E_coulomb) ;
            Atoms[i].Pres_Work+= 0.5*r_ij*Force;
			Atoms[j].Pres_Work+= 0.5*r_ij*Force;
				
			if ((Atoms[i].R_Vec[0]>X_L_Measure)&&(Atoms[i].R_Vec[0]<X_H_Measure)) {
				if ((Atoms[i].R_Vec[1]>Y_L_Measure)&&(Atoms[i].R_Vec[1]<Y_H_Measure)) {
					if ((Atoms[i].R_Vec[2]>Z_L_Measure)&&(Atoms[i].R_Vec[2]<Z_H_Measure)) {
						W_SR+=r_ij2*Force;							
					}
				}
			}
		}
	
	}

	for (i=0;i<Num_Of_Atoms;i++) { 
		Pot_tot += Atoms[i].Pot ;
	}
}
/***********************************************************************/
/***********************************************************************/
void LJ_COUL_of_Atom(int i) 
{
	int	j,jj;
	double x_ij,y_ij,z_ij,r_ij2,r_ij2_inv,r_ij6_inv,r_ij12_inv,Force,r_ij,x_ij2,y_ij2,z_ij2;
	double E_coulomb,E_LJ=0, Charge1, Charge2, rinv, forcecoul, r_inv_ij2, LJ_Force, COUL_Force, RR_Cut_Off_R;
	int numofNeighbors;
	double sigma,epsilon;
	numofNeighbors=Atoms[i].NeighborList[0];
	E_coulomb = 0.0;
	Atoms[i].Pot=0;
	for (jj=1;jj<numofNeighbors+1;jj++) {
		j=Atoms[i].NeighborList[jj];
		z_ij=(Atoms[i].R_Vec[2]-Atoms[j].R_Vec[2]);
		y_ij=(Atoms[i].R_Vec[1]-Atoms[j].R_Vec[1]);
		x_ij=(Atoms[i].R_Vec[0]-Atoms[j].R_Vec[0]);
		if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
			if ((x_ij)>X_H_L/2)  x_ij=x_ij-X_H_L;
			else if ((x_ij)<-X_H_L/2) x_ij=x_ij+X_H_L;
		}
		if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
			if ((y_ij)>Y_H_L/2) y_ij=y_ij-Y_H_L;
			else if ((y_ij)<-Y_H_L/2) y_ij=y_ij+Y_H_L;
		}
		if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
			if ((z_ij)>Z_H_L/2) z_ij=z_ij-Z_H_L;
			else if ((z_ij)<-Z_H_L/2) z_ij=z_ij+Z_H_L;
		}
		x_ij2 = x_ij*x_ij;
		y_ij2 = y_ij*y_ij;
		z_ij2 = z_ij*z_ij;
		r_ij2 = (x_ij2 + y_ij2 + z_ij2);
		r_ij = sqrt(r_ij2);
		sigma = Sigma_MAT[Atoms[i].tp-1][Atoms[j].tp-1];
		R_Cut_Off = R_Cut_Off_R*sigma;
		Charge1 = q[Atoms[i].tp-1];
		Charge2 = q[Atoms[j].tp-1];
		RR_Cut_Off_R = R_Cut_Off*R_Cut_Off;
		COUL_Force = 0.0;
		E_coulomb=0;
	    if (fabs(Charge1) >0.000001) {
			if (fabs(Charge2) >0.000001) {
				if ((r_ij<R_Cut_Coul)) {	
					r_inv_ij2 = (1.0)/r_ij2;
					rinv = sqrt(r_inv_ij2);
					forcecoul = -dielectric_constant * Charge1 * Charge2 * rinv;
					COUL_Force = (coulomb_factor*r_inv_ij2)*(-forcecoul);
					E_coulomb = coulomb_factor * (forcecoul);
				}
			}
		}
		LJ_Force = 0.0;
		E_LJ=0;
		if ((r_ij<R_Cut_Off)) {
			r_ij2_inv = (sigma*sigma)/r_ij2;
			r_ij6_inv = r_ij2_inv*r_ij2_inv*r_ij2_inv;
			r_ij12_inv = r_ij6_inv*r_ij6_inv;
			r_ij2_inv = r_ij2_inv/(sigma*sigma);
			epsilon = Epsilon_MAT[Atoms[i].tp-1][Atoms[j].tp-1];
			LJ_Force = 48*epsilon*r_ij2_inv*(r_ij12_inv-0.5*r_ij6_inv);
			E_LJ=4*epsilon*(r_ij12_inv-r_ij6_inv);
		}
		Force = LJ_Force + COUL_Force;
		Atoms[i].force[0] += x_ij*Force;
		Atoms[i].force[1] += y_ij*Force;
		Atoms[i].force[2] += z_ij*Force;
		Atoms[j].force[0] -= x_ij*Force;
		Atoms[j].force[1] -= y_ij*Force;
		Atoms[j].force[2] -= z_ij*Force;
		Atoms[i].Pot +=0.5*(E_LJ+E_coulomb) ;
		Atoms[j].Pot +=0.5*(E_LJ+E_coulomb) ;
        Atoms[i].Pres_Work+= 0.5*r_ij*Force;
		Atoms[j].Pres_Work+= 0.5*r_ij*Force;
				
		if ((Atoms[i].R_Vec[0]>X_L_Measure)&&(Atoms[i].R_Vec[0]<X_H_Measure)) {
			if ((Atoms[i].R_Vec[1]>Y_L_Measure)&&(Atoms[i].R_Vec[1]<Y_H_Measure)) {
				if ((Atoms[i].R_Vec[2]>Z_L_Measure)&&(Atoms[i].R_Vec[2]<Z_H_Measure)) {
					W_SR+=r_ij2*Force;							
				}
			}
		}
	}

	for (i=0;i<Num_Of_Atoms;i++) { 
		Pot_tot += Atoms[i].Pot ;
	}
}
/***********************************************************************/
/***********************************************************************/
void Bond_Force()
{
	Harmonic_Bond();
}
/***********************************************************************/
/***********************************************************************/
void Harmonic_Bond()
{
  int i1,i2,n,tp;
  double fbond,E_Bound=0;
  double x_ij,y_ij,z_ij,r_ij,dr,r01,rk;
  
	  for (n = 0; n < Num_Of_Bonds; n++) {
		 i1 = bondlist[n][0];
		 i2 = bondlist[n][1];
		 tp = bondlist[n][2];
		 x_ij = Atoms[i1].R_Vec[0] - Atoms[i2].R_Vec[0];
		 y_ij = Atoms[i1].R_Vec[1] - Atoms[i2].R_Vec[1];
		 z_ij = Atoms[i1].R_Vec[2] - Atoms[i2].R_Vec[2];
   
		if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
			if ((x_ij)>(X_High-X_Low)/2) x_ij=x_ij-(X_High-X_Low);
			else if ((x_ij)<(X_Low-X_High)/2) x_ij=x_ij+(X_High-X_Low);
		}
		if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
			if ((y_ij)>(Y_High-Y_Low)/2) y_ij=y_ij-(Y_High-Y_Low);
			else if ((y_ij)<(Y_Low-Y_High)/2) y_ij=y_ij+(Y_High-Y_Low);
		}
		
		if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
			if ((z_ij)>(Z_High-Z_Low)/2) z_ij=z_ij-(Z_High-Z_Low);
			else if ((z_ij)<(Z_Low-Z_High)/2) z_ij=z_ij+(Z_High-Z_Low);
		}
		
		r_ij=sqrt(x_ij*x_ij + y_ij*y_ij + z_ij*z_ij);
		r01=Harmonic_Bond_Parameter[tp-1].r0;
		dr=r_ij-r01;
		rk=dr*Harmonic_Bond_Parameter[tp-1].K;
		if (r_ij > 0.0) fbond = -rk/r_ij; 
		else			fbond = 0.0;

		Atoms[i1].force[0]+= x_ij*fbond;
		Atoms[i1].force[1]+= y_ij*fbond;
		Atoms[i1].force[2]+= z_ij*fbond;
		Atoms[i2].force[0]-= x_ij*fbond;
		Atoms[i2].force[1]-= y_ij*fbond;
		Atoms[i2].force[2]-= z_ij*fbond;

		E_Bound+=rk*dr;
	  }
	  Pot_tot+=0.5*E_Bound;
}
/***********************************************************************/
/***********************************************************************/
void Const_Temp_and_Flux_Wall_Force()
{
  int i,j,J;
  double fT,ET=0,U,V,W;
  double x_ij,y_ij,z_ij,r_ij,rk;
  /*********************** CONSTANT TEMPERATURE ************************/
  for (i=0; i< Num_Of_Cons_Temp_Sets; i++) {
	U=cons_Temp_Wall[i].vx;
	V=cons_Temp_Wall[i].vy;
	W=cons_Temp_Wall[i].vz;
	for (j =0; j < Cons_Temp_Atoms[i][0]; j++) {
		J=Cons_Temp_Atoms[i][j+1];
		 x_ij = Atoms[J].R_Vec[0] - X_e_T_Wall[i][j];
	     y_ij = Atoms[J].R_Vec[1] - Y_e_T_Wall[i][j];
		 z_ij = Atoms[J].R_Vec[2] - Z_e_T_Wall[i][j];
		if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
			if ((x_ij)>(X_H_L)/2) x_ij=x_ij-(X_H_L);
			else if ((x_ij)<(-X_H_L)/2) x_ij=x_ij+(X_H_L);
		}
		if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
			if ((y_ij)>(Y_H_L)/2) y_ij=y_ij-(Y_H_L);
			else if ((y_ij)<(-Y_H_L)/2) y_ij=y_ij+(Y_H_L);
		}
		if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
			if ((z_ij)>(Z_H_L)/2) z_ij=z_ij-(Z_H_L);
			else if ((z_ij)<(-Z_H_L)/2) z_ij=z_ij+(Z_H_L);
		}
		r_ij=sqrt(x_ij*x_ij+y_ij*y_ij+z_ij*z_ij);
		rk=r_ij*cons_Temp_Wall[i].K;
		if (r_ij > 0.0) fT = -rk/r_ij; 
		else fT = 0.0;
		Atoms[J].force[0]+= x_ij*fT;
		Atoms[J].force[1]+= y_ij*fT;
		Atoms[J].force[2]+= z_ij*fT;
        X_e_T_Wall[i][j]+=U*DeltaT;
        Y_e_T_Wall[i][j]+=V*DeltaT;
		Z_e_T_Wall[i][j]+=W*DeltaT;
		
		if (X_e_T_Wall[i][j]<X_Low) {
			X_e_T_Wall[i][j]=X_H_L+X_e_T_Wall[i][j];
        } else	if (X_e_T_Wall[i][j]>X_High) {
	   		X_e_T_Wall[i][j]=-X_H_L+X_e_T_Wall[i][j];
	    }
		if (Y_e_T_Wall[i][j]<Y_Low) {
	   		Y_e_T_Wall[i][j]=Y_H_L+Y_e_T_Wall[i][j];
	    } else	if (Y_e_T_Wall[i][j]>Y_High) {
	   		Y_e_T_Wall[i][j]=-Y_H_L+Y_e_T_Wall[i][j];
	    }
		if (Z_e_T_Wall[i][j]<Z_Low) {
	   		Z_e_T_Wall[i][j]=Z_H_L+Z_e_T_Wall[i][j];
	    } else	if (Z_e_T_Wall[i][j]>Z_High) {
	   		Z_e_T_Wall[i][j]=-Z_H_L+Z_e_T_Wall[i][j];
	    }
		ET+=rk*r_ij;
	}
  }
  /*********************** CONSTANT FLUX ************************/
  for (i=0; i< Num_Of_Const_Flux_Sets; i++) {
  	  U=cons_Flux_Wall[i].vx;
	  V=cons_Flux_Wall[i].vy;
	  W=cons_Flux_Wall[i].vz;
	  for (j =0; j < Const_Flux_Atoms[i][0]; j++) {	
		J=Const_Flux_Atoms[i][j+1];
		 x_ij = Atoms[J].R_Vec[0] - X_e_F_Wall[i][j];
	     y_ij = Atoms[J].R_Vec[1] - Y_e_F_Wall[i][j];
		 z_ij = Atoms[J].R_Vec[2] - Z_e_F_Wall[i][j];
		if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
 			if ((x_ij)>(X_H_L)/2) x_ij=x_ij-(X_H_L);
			else if ((x_ij)<(-X_H_L)/2) x_ij=x_ij+(X_H_L);
		}
		if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
			if ((y_ij)>(Y_H_L)/2) y_ij=y_ij-(Y_H_L);
			else if ((y_ij)<(-Y_H_L)/2) y_ij=y_ij+(Y_H_L);
		}
		if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
			if ((z_ij)>(Z_H_L)/2) z_ij=z_ij-(Z_H_L);
			else if ((z_ij)<(-Z_H_L)/2) z_ij=z_ij+(Z_H_L);
		}
		r_ij=sqrt(x_ij*x_ij+y_ij*y_ij+z_ij*z_ij);
		rk=r_ij*cons_Flux_Wall[i].K;
		if (r_ij > 0.0) fT = -rk/r_ij; 
		else fT = 0.0;
		Atoms[J].force[0]+= x_ij*fT;
		Atoms[J].force[1]+= y_ij*fT;
		Atoms[J].force[2]+= z_ij*fT;
        X_e_F_Wall[i][j]+=U*DeltaT;
        Y_e_F_Wall[i][j]+=V*DeltaT;
		Z_e_F_Wall[i][j]+=W*DeltaT;
		if (X_e_F_Wall[i][j]<X_Low) {
    		X_e_F_Wall[i][j]=X_H_L+X_e_F_Wall[i][j];
         } else	if (X_e_F_Wall[i][j]>X_High) {
     		X_e_F_Wall[i][j]=-X_H_L+X_e_F_Wall[i][j];
         }
		if (Y_e_F_Wall[i][j]<Y_Low) {
    		Y_e_F_Wall[i][j]=Y_H_L+Y_e_F_Wall[i][j];
        } else	if (Y_e_F_Wall[i][j]>Y_High) {
     		Y_e_F_Wall[i][j]=-Y_H_L+Y_e_F_Wall[i][j];
        }
		if (Z_e_F_Wall[i][j]<Z_Low) {
    		Z_e_F_Wall[i][j]=Z_H_L+Z_e_F_Wall[i][j];
        } else	if (Z_e_F_Wall[i][j]>Z_High) {
     		Z_e_F_Wall[i][j]=-Z_H_L+Z_e_F_Wall[i][j];
        }
		ET+=rk*r_ij;
	}
 }
	Pot_tot+=0.5*ET;
}
/***********************************************************************/
/***********************************************************************/
void Angle_Force()
{
	Harmonic_Angle();
}
/***********************************************************************/
/***********************************************************************/
void Harmonic_Angle()
{
	int i1,i2,i3,n,tp;
	double x_ij,y_ij,z_ij,x_kj,y_kj,z_kj;
	double E_angle=0,f1[3],f3[3];
	double tk,dtheta,r2_ij;
	double r_ij,r_kj,c,s,a,a11,a12,a22,r2_kj;
	double pi_N=3.1415927;
	double E_Angle=0;
	for (n = 0; n < Num_Of_Angles; n++) {
		i1 = anglelist[n][0];
		i2 = anglelist[n][1];
		i3 = anglelist[n][2];
		tp = anglelist[n][3];
		/*************** 1st BOND  *******************/
		x_ij = Atoms[i1].R_Vec[0] - Atoms[i2].R_Vec[0];
		y_ij = Atoms[i1].R_Vec[1] - Atoms[i2].R_Vec[1];
		z_ij = Atoms[i1].R_Vec[2] - Atoms[i2].R_Vec[2];
		if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
			if ((x_ij)>(X_High-X_Low)/2) x_ij=x_ij-(X_High-X_Low);
			else if ((x_ij)<(X_Low-X_High)/2) x_ij=x_ij+(X_High-X_Low);
		}
		if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
			if ((y_ij)>(Y_High-Y_Low)/2) y_ij=y_ij-(Y_High-Y_Low);
			else if ((y_ij)<(Y_Low-Y_High)/2) y_ij=y_ij+(Y_High-Y_Low);
		}
	
		if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
			if ((z_ij)>(Z_High-Z_Low)/2) z_ij=z_ij-(Z_High-Z_Low);
			else if ((z_ij)<(Z_Low-Z_High)/2) z_ij=z_ij+(Z_High-Z_Low);
		}
	
		r2_ij=(x_ij*x_ij+y_ij*y_ij+z_ij*z_ij);
		r_ij=sqrt(r2_ij);
		/*************** 2nd BOND  *******************/
		x_kj = Atoms[i3].R_Vec[0] - Atoms[i2].R_Vec[0];
		y_kj = Atoms[i3].R_Vec[1] - Atoms[i2].R_Vec[1];
		z_kj = Atoms[i3].R_Vec[2] - Atoms[i2].R_Vec[2];
    
		if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
			if ((x_kj)>(X_High-X_Low)/2) x_kj=x_kj-(X_High-X_Low);
			else if ((x_kj)<(X_Low-X_High)/2) x_kj=x_kj+(X_High-X_Low);
		}
		if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
			if ((y_kj)>(Y_High-Y_Low)/2) y_kj=y_kj-(Y_High-Y_Low);
			else if ((y_kj)<(Y_Low-Y_High)/2) y_kj=y_kj+(Y_High-Y_Low);
		}
	
		if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
			if ((z_kj)>(Z_High-Z_Low)/2) z_kj=z_kj-(Z_High-Z_Low);
			else if ((z_kj)<(Z_Low-Z_High)/2) z_kj=z_kj+(Z_High-Z_Low);
		}

		r2_kj=(x_kj*x_kj+y_kj*y_kj+z_kj*z_kj);
		r_kj=sqrt(r2_kj);
		c = x_ij*x_kj + y_ij*y_kj + z_ij*z_kj;
		c /= r_ij*r_kj;
		if (c > 1.0) c = 1.0;
		if (c < -1.0) c = -1.0;
	    s = sqrt(1.0 - c*c);
		s = 1.0/s;

    /***************** FORCE & ENERGY ******************/
		dtheta = acos(c) - Harmonic_Angle_Parameter[tp-1].r0;
        tk = Harmonic_Angle_Parameter[tp-1].K * dtheta;
        a = -tk * s;
        a11 = a*c / r2_ij;
        a12 = -a / (r_ij*r_kj);
        a22 = a*c / r2_kj;
	    E_Angle+= tk*dtheta;
		f1[0] = a11*x_ij + a12*x_kj;
		f1[1] = a11*y_ij + a12*y_kj;
		f1[2] = a11*z_ij + a12*z_kj;
		f3[0] = a22*x_kj + a12*x_ij;
		f3[1] = a22*y_kj + a12*y_ij;
		f3[2] = a22*z_kj + a12*z_ij;
		Atoms[i1].force[0] += f1[0];
		Atoms[i1].force[1] += f1[1];
		Atoms[i1].force[2] += f1[2];
   		Atoms[i2].force[0] -= f1[0] + f3[0];
		Atoms[i2].force[1] -= f1[1] + f3[1];
		Atoms[i2].force[2] -= f1[2] + f3[2];
		Atoms[i3].force[0] += f3[0];
		Atoms[i3].force[1] += f3[1];
		Atoms[i3].force[2] += f3[2];
    }
	Pot_tot+=0.5*E_Angle;
}
/***********************************************************************/
/***********************************************************************/
void Dihedral_Force() 
{
	Harmonic_Dihedral();
}
/***********************************************************************/
/***********************************************************************/
void Harmonic_Dihedral()
{
  int i1,i2,i3,i4,i,m,n,tp;
  double x_ij,y_ij,z_ij,x_kj,y_kj,z_kj,x_lj,y_lj,z_lj,xm_kj,ym_kj,zm_kj;
  double E_dihedral,f1[3],f2[3],f3[3],f4[3];
  double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,rginv,ra2inv,rb2inv,rabinv;
  double df,df1,ddf1,fg,hg,fga,hgb,gaa,gbb;
  double dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx,dthy,dthz;  
  double c,s,p,sx2,sy2,sz2;

  E_dihedral = 0.0;

  // require sign = +/- 1 for backwards compatibility
  // arbitrary phase angle shift could be allowed, but would break
  //   backwards compatibility and is probably not needed
                       
  int count = 0;
  for (i = 0; i <Num_Of_Dihedral_Types ; i++) {

	  if (Harmonic_Dihedral_Parameter[i].sign != -1 && Harmonic_Dihedral_Parameter[i].sign != 1)	cout << "Incorrect sign arg for dihedral coefficients" << endl;
	  if (Harmonic_Dihedral_Parameter[i].multiplicity < 0)	cout << "Incorrect multiplicity arg for dihedral coefficients" << endl;
     
	  if (Harmonic_Dihedral_Parameter[i].sign == 1) {
		Harmonic_Dihedral_Parameter[i].cos_shift = 1;
		Harmonic_Dihedral_Parameter[i].sin_shift = 0;
	  } 
	  else 
	  {
      Harmonic_Dihedral_Parameter[i].cos_shift = -1;
      Harmonic_Dihedral_Parameter[i].sin_shift = 0;
	  }
      count++;
  } 

  for (n = 0; n < Num_Of_Dihedrals; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    tp = dihedrallist[n][4];

    /*************** 1st BOND  *******************/

    x_ij = Atoms[i1].R_Vec[0] - Atoms[i2].R_Vec[0];
    y_ij = Atoms[i1].R_Vec[1] - Atoms[i2].R_Vec[1];
    z_ij = Atoms[i1].R_Vec[2] - Atoms[i2].R_Vec[2];

    if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
    	if ((x_ij)>(X_High-X_Low)/2) x_ij=x_ij-(X_High-X_Low);
		else if ((x_ij)<(X_Low-X_High)/2) x_ij=x_ij+(X_High-X_Low);
    }
    if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
		if ((y_ij)>(Y_High-Y_Low)/2) y_ij=y_ij-(Y_High-Y_Low);
		else if ((y_ij)<(Y_Low-Y_High)/2) y_ij=y_ij+(Y_High-Y_Low);
    }
    if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
		if ((z_ij)>(Z_High-Z_Low)/2) z_ij=z_ij-(Z_High-Z_Low);
		else if ((z_ij)<(Z_Low-Z_High)/2) z_ij=z_ij+(Z_High-Z_Low);
    }
    /*************** 2nd BOND  *******************/
    x_kj = Atoms[i2].R_Vec[0] - Atoms[i2].R_Vec[0];
    y_kj = Atoms[i3].R_Vec[1] - Atoms[i2].R_Vec[1];
    z_kj = Atoms[i3].R_Vec[2] - Atoms[i2].R_Vec[2];

    if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
		if ((x_kj)>(X_High-X_Low)/2) x_kj=x_kj-(X_High-X_Low);
		else if ((x_kj)<(X_Low-X_High)/2) x_kj=x_kj+(X_High-X_Low);
    }
    if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
        if ((y_kj)>(Y_High-Y_Low)/2) y_kj=y_kj-(Y_High-Y_Low);
		else if ((y_kj)<(Y_Low-Y_High)/2) y_kj=y_kj+(Y_High-Y_Low);
    }
    if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
		if ((z_kj)>(Z_High-Z_Low)/2) z_kj=z_kj-(Z_High-Z_Low);
		else if ((z_kj)<(Z_Low-Z_High)/2) z_kj=z_kj+(Z_High-Z_Low);
    }

    xm_kj = -x_kj;
    ym_kj = -y_kj;
    zm_kj = -z_kj;

    if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
		if ((xm_kj)>(X_High-X_Low)/2) xm_kj=xm_kj-(X_High-X_Low);
		else if ((xm_kj)<(X_Low-X_High)/2) xm_kj=xm_kj+(X_High-X_Low);
    }
    if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
        if ((ym_kj)>(Y_High-Y_Low)/2) ym_kj=ym_kj-(Y_High-Y_Low);
		else if ((ym_kj)<(Y_Low-Y_High)/2) ym_kj=ym_kj+(Y_High-Y_Low);
    }
   if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
		if ((zm_kj)>(Z_High-Z_Low)/2) zm_kj=zm_kj-(Z_High-Z_Low);
		else if ((zm_kj)<(Z_Low-Z_High)/2) zm_kj=zm_kj+(Z_High-Z_Low);
    }

    /*************** 3rd BOND  *******************/
    x_lj = Atoms[i4].R_Vec[0] - Atoms[i3].R_Vec[0];
    y_lj = Atoms[i4].R_Vec[1] - Atoms[i3].R_Vec[1];
    z_lj = Atoms[i4].R_Vec[2] - Atoms[i3].R_Vec[2];
 
    if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
		if ((x_lj)>(X_High-X_Low)/2) x_lj=x_lj-(X_High-X_Low);
		else if ((x_lj)<(X_Low-X_High)/2) x_lj=x_lj+(X_High-X_Low);
    }
    if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
        if ((y_lj)>(Y_High-Y_Low)/2) y_lj=y_lj-(Y_High-Y_Low);
		else if ((y_lj)<(Y_Low-Y_High)/2) y_lj=y_lj+(Y_High-Y_Low);
    }
   if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
		if ((z_lj)>(Z_High-Z_Low)/2) z_lj=z_lj-(Z_High-Z_Low);
		else if ((z_lj)<(Z_Low-Z_High)/2) z_lj=z_lj+(Z_High-Z_Low);
    }

    // SIN And COSINE
    ax = y_ij*zm_kj - z_ij*ym_kj;
    ay = z_ij*xm_kj - x_ij*zm_kj;
    az = x_ij*ym_kj - y_ij*xm_kj;
    bx = y_lj*zm_kj - z_lj*ym_kj;
    by = z_lj*xm_kj - x_lj*zm_kj;
    bz = x_lj*ym_kj - y_lj*xm_kj;

    rasq = ax*ax + ay*ay + az*az;
    rbsq = bx*bx + by*by + bz*bz;
    rgsq = xm_kj*xm_kj + ym_kj*ym_kj + zm_kj*zm_kj;
    rg = sqrt(rgsq);
    
    rginv = ra2inv = rb2inv = 0.0;
    if (rg > 0) rginv = 1.0/rg;
    if (rasq > 0) ra2inv = 1.0/rasq;
    if (rbsq > 0) rb2inv = 1.0/rbsq;
    rabinv = sqrt(ra2inv*rb2inv);

    c = (ax*bx + ay*by + az*bz)*rabinv;
    s = rg*rabinv*(ax*x_lj + ay*y_lj + az*z_lj);
    
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
         
    m = Harmonic_Dihedral_Parameter[tp-1].multiplicity;
    p = 1.0;
    df1 = 0.0;
    
    for (i = 0; i < m; i++) {
      ddf1 = p*c - df1*s;
      df1 = p*s + df1*c;
      p = ddf1;
    }

    p = p * Harmonic_Dihedral_Parameter[tp-1].cos_shift + df1 * Harmonic_Dihedral_Parameter[tp-1].sin_shift;
    df1 = df1 * Harmonic_Dihedral_Parameter[tp-1].cos_shift - ddf1 * Harmonic_Dihedral_Parameter[tp-1].sin_shift;
    df1 *= -m;
    p += 1.0;
 
    if (m == 0) {
      p = 1.0 + Harmonic_Dihedral_Parameter[tp-1].cos_shift;
      df1 = 0.0;
    }

    E_dihedral += Harmonic_Dihedral_Parameter[tp-1].K * p; 
       
    fg = x_ij*xm_kj + y_ij*ym_kj + z_ij*zm_kj;
    hg = x_lj*xm_kj + y_lj*ym_kj + z_lj*zm_kj;
    fga = fg*ra2inv*rginv;
    hgb = hg*rb2inv*rginv;
    gaa = -ra2inv*rg;
    gbb = rb2inv*rg;
    
    dtfx = gaa*ax;
    dtfy = gaa*ay;
    dtfz = gaa*az;
    dtgx = fga*ax - hgb*bx;
    dtgy = fga*ay - hgb*by;
    dtgz = fga*az - hgb*bz;
    dthx = gbb*bx;
    dthy = gbb*by;
    dthz = gbb*bz;
    
    df = -Harmonic_Dihedral_Parameter[tp-1].K * df1;
    
    sx2 = df*dtgx;
    sy2 = df*dtgy;
    sz2 = df*dtgz;

    f1[0] = df*dtfx;
    f1[1] = df*dtfy;
    f1[2] = df*dtfz;

    f2[0] = sx2 - f1[0];
    f2[1] = sy2 - f1[1];
    f2[2] = sz2 - f1[2];

    f4[0] = df*dthx;
    f4[1] = df*dthy;
    f4[2] = df*dthz;

    f3[0] = -sx2 - f4[0];
    f3[1] = -sy2 - f4[1];
    f3[2] = -sz2 - f4[2];
    
    Atoms[i1].force[0] += f1[0];
    Atoms[i1].force[1] += f1[1];
    Atoms[i1].force[2] += f1[2];
    Atoms[i2].force[0] += f2[0];
    Atoms[i2].force[1] += f2[1];
    Atoms[i2].force[2] += f2[2];
   
    Atoms[i3].force[0] += f3[0];
    Atoms[i3].force[1] += f3[1];
    Atoms[i3].force[2] += f3[2];
    
    Atoms[i4].force[0] += f4[0];
    Atoms[i4].force[1] += f4[1];
    Atoms[i4].force[2] += f4[2];
  }
  Pot_tot += E_dihedral;
}
/***********************************************************************/
/***********************************************************************/
void Integrate_V()
{
	int		i;
	double dt2=DeltaT*DeltaT/2;
	double dt=DeltaT/2;
	double delx,dely,delz,delr,delr_max,a_mass;
	double fact;
	xip=new double[Num_Of_Cons_Temp_Sets];
	delr_max=0;
	/**************************** THERMOSTAT ***************************/
	if (Thermo_tp=='Y') {
		fact=(1-DeltaT*Xi/2);
		for (i=0;i<Num_Of_Atoms;i++) {
			Atoms[i].V_Vec[0]*=fact;
			Atoms[i].V_Vec[1]*=fact;
			Atoms[i].V_Vec[2]*=fact;
		}
	}
	/******************************************************************/
	for (i=0;i<Num_Of_Atoms;i++) {
		delx=Atoms[i].V_Vec[0]*DeltaT+Atoms[i].A_Vec[0]*(dt2);
		dely=Atoms[i].V_Vec[1]*DeltaT+Atoms[i].A_Vec[1]*(dt2);
		delz=Atoms[i].V_Vec[2]*DeltaT+Atoms[i].A_Vec[2]*(dt2);
		Atoms[i].R_Vec[0]+=delx;
		Atoms[i].R_Vec[1]+=dely;
		Atoms[i].R_Vec[2]+=delz;
		delr=sqrt(delx*delx+dely*dely+delz*delz);
		if (delr>delr_max) delr_max=delr;
	}
	
	for (i=0;i<Num_Of_Atoms;i++) {
		Atoms[i].V_Vec[0]+=Atoms[i].A_Vec[0]*dt;
		Atoms[i].V_Vec[1]+=Atoms[i].A_Vec[1]*dt;
		Atoms[i].V_Vec[2]+=Atoms[i].A_Vec[2]*dt;
	}
	if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P'))	{
		for (i=0;i<Num_Of_Atoms;i++) 		X_Periodic_Boundary(i) ;
	}
	if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
		for (i=0;i<Num_Of_Atoms;i++) Y_Periodic_Boundary(i) ;
	}
	if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
		for (i=0;i<Num_Of_Atoms;i++)   Z_Periodic_Boundary(i) ;
	}
	rebuild_neighb_char+=delr_max;
	Compute_Force();
	for (i=0;i<Num_Of_Atoms;i++) {
		a_mass=AtomMass[Atoms[i].tp-1];
		Atoms[i].A_Vec[0]=Atoms[i].force[0]/a_mass;
		Atoms[i].A_Vec[1]=Atoms[i].force[1]/a_mass;
		Atoms[i].A_Vec[2]=Atoms[i].force[2]/a_mass;
	}
	
	for (i=0;i<Num_Of_Atoms;i++) {
		Atoms[i].V_Vec[0]+=Atoms[i].A_Vec[0]*dt;
		Atoms[i].V_Vec[1]+=Atoms[i].A_Vec[1]*dt;
		Atoms[i].V_Vec[2]+=Atoms[i].A_Vec[2]*dt;
	}
	/////////////// 
	if (Thermo_tp=='Y') {
		fact=(1+DeltaT*Xi/2);
		for (i=0;i<Num_Of_Atoms;i++) {
			Atoms[i].V_Vec[0]=Atoms[i].V_Vec[0]/fact;
			Atoms[i].V_Vec[1]=Atoms[i].V_Vec[1]/fact;
			Atoms[i].V_Vec[2]=Atoms[i].V_Vec[2]/fact;
		}
		Xip=Xi;
    	Calc_Properties_Mid() ;
	    fact=(1+DeltaT*Xip/2)/(1+DeltaT*Xi/2);	
		for (i=0;i<Num_Of_Atoms;i++) {
			Atoms[i].V_Vec[0]=Atoms[i].V_Vec[0]*fact;
			Atoms[i].V_Vec[1]=Atoms[i].V_Vec[1]*fact;
			Atoms[i].V_Vec[2]=Atoms[i].V_Vec[2]*fact;
		}
		Xip=Xi;
		Calc_Properties_Mid() ;
		fact=(1+DeltaT*Xip/2)/(1+DeltaT*Xi/2);
		for (i=0;i<Num_Of_Atoms;i++) {
			Atoms[i].V_Vec[0]=Atoms[i].V_Vec[0]*fact;
			Atoms[i].V_Vec[1]=Atoms[i].V_Vec[1]*fact;
			Atoms[i].V_Vec[2]=Atoms[i].V_Vec[2]*fact;
		}
		Xip=Xi;
		Calc_Properties_Mid() ;
		fact=(1+DeltaT*Xip/2)/(1+DeltaT*Xi/2);
		
		for (i=0;i<Num_Of_Atoms;i++) {
			Atoms[i].V_Vec[0]=Atoms[i].V_Vec[0]*fact;
			Atoms[i].V_Vec[1]=Atoms[i].V_Vec[1]*fact;
			Atoms[i].V_Vec[2]=Atoms[i].V_Vec[2]*fact;
		}
	}
	Fix_Velocity() ;
	if (Wall_Thermo_Method=='S') Wall_Temp_Scaling();
	delete [] xip;
	xip=NULL;
	Delete_Bound_OK=0;
	if ((Boundary_Type.x=='f') |(Boundary_Type.x=='F'))	{
		for (i=0;i<Num_Of_Atoms;i++) 		X_Free_Boundary(i) ;
	}
	if ((Boundary_Type.y=='f') |(Boundary_Type.y=='F')) {
		for (i=0;i<Num_Of_Atoms;i++) Y_Free_Boundary(i) ;
	}
	if ((Boundary_Type.z=='f') |(Boundary_Type.z=='F')) {
		for (i=0;i<Num_Of_Atoms;i++)   Z_Free_Boundary(i) ;
	}
	if (Delete_Bound_OK==1) {
		Delete_Atoms();
		Ref_Del_List();
	}
}
/***********************************************************************/
/***********************************************************************/
void	Calc_Properties_Mid() 
{
	int	i,i1;
	double	Kin_tot_Measure=0;
	double	M_sys=0,a_mass;
	N_Measure=0;
	v_mean[0]=0;
	v_mean[1]=0;
	v_mean[2]=0;
	for (i=0;i<Num_Of_Atoms;i++) {
		Atoms[i].Inc_vec=0;
		if ((Atoms[i].R_Vec[0]>X_L_Measure)&&(Atoms[i].R_Vec[0]<X_H_Measure)) {
			if ((Atoms[i].R_Vec[1]>Y_L_Measure)&&(Atoms[i].R_Vec[1]<Y_H_Measure)) {
				if ((Atoms[i].R_Vec[2]>Z_L_Measure)&&(Atoms[i].R_Vec[2]<Z_H_Measure)) {
				//	Inc=1;
					Atoms[i].Inc_vec=1;
					for (i1=0;i1<Num_Of_Vel_Sets;i1++) {
						if ((Atoms[i].tp==Vel_Set[i1].tp)) Atoms[i].Inc_vec=0;
					}
					for (i1=0;i1<Num_Of_Cons_Temp_Sets;i1++) {
						if (Atoms[i].tp==cons_Temp_Wall[i1].tp) Atoms[i].Inc_vec=0;
					}
					if (Atoms[i].Inc_vec==1) {
						a_mass=AtomMass[Atoms[i].tp-1];
						v_mean[0]+=(a_mass*Atoms[i].V_Vec[0]);
						v_mean[1]+=(a_mass*Atoms[i].V_Vec[1]);
						v_mean[2]+=(a_mass*Atoms[i].V_Vec[2]);
						N_Measure++;
						M_sys+=a_mass;
					}
				}
			}
		}
	}
	v_mean[0]=v_mean[0]/M_sys;
	v_mean[1]=v_mean[1]/M_sys;
	v_mean[2]=v_mean[2]/M_sys;
	for (i=0;i<Num_Of_Atoms;i++) {
		if (Atoms[i].Inc_vec==1) {
			Kin_tot_Measure+=AtomMass[Atoms[i].tp-1]*((Atoms[i].V_Vec[0]-v_mean[0])*(Atoms[i].V_Vec[0]-v_mean[0])+(Atoms[i].V_Vec[1]-v_mean[1])*(Atoms[i].V_Vec[1]-v_mean[1])+(Atoms[i].V_Vec[2]-v_mean[2])*(Atoms[i].V_Vec[2]-v_mean[2]));
		}
	}
	Kin_tot_Measure=Kin_tot_Measure/2;
	if (Thermo_tp=='Y') {
		Xi+=DeltaT*(Kin_tot_Measure-1.5*N_Measure*Boltzmann_const*Target_Temp)/T_Coupling;
	}
}
/***********************************************************************/
/***********************************************************************/
void	Calc_Properties() 
{
	int	i;
	double	Volume_Measure;
	double	Kin_tot_Measure=0;
	double	M_sys=0;
	Kin_tot=0;
	for (i=0;i<Num_Of_Atoms;i++) {
		if (Atoms[i].Inc_vec==1) {
			Kin_tot_Measure+=AtomMass[Atoms[i].tp-1]*((Atoms[i].V_Vec[0]-v_mean[0])*(Atoms[i].V_Vec[0]-v_mean[0])+(Atoms[i].V_Vec[1]-v_mean[1])*(Atoms[i].V_Vec[1]-v_mean[1])+(Atoms[i].V_Vec[2]-v_mean[2])*(Atoms[i].V_Vec[2]-v_mean[2]));
		}
		Kin_tot+=AtomMass[Atoms[i].tp-1]*((Atoms[i].V_Vec[0])*(Atoms[i].V_Vec[0])+(Atoms[i].V_Vec[1])*(Atoms[i].V_Vec[1])+(Atoms[i].V_Vec[2])*(Atoms[i].V_Vec[2]));
	}
	Kin_tot_Measure=Kin_tot_Measure/2;
	Kin_tot=Kin_tot/2;
	V_mean=sqrt(v_mean[0]*v_mean[0]+v_mean[1]*v_mean[1]+v_mean[2]*v_mean[2]);
	Volume_Measure=(X_H_Measure-X_L_Measure)*(Y_H_Measure-Y_L_Measure)*(Z_H_Measure-Z_L_Measure);
	Sys_Temperature=Kin_tot_Measure/(1.5*N_Measure*Boltzmann_const);
	Sys_Pressure=(Kin_tot_Measure+W_SR/2.0)/(1.5*Volume_Measure);
}
/***********************************************************************/
/***********************************************************************/
void	Show_Step_Prop()
{
	if (StepNum%Prop_Write_Step==0) {
		pot_file << StepNum << "      "  << Time << "  " << Pot_tot/Num_Of_Atoms << "  " << Kin_tot/Num_Of_Atoms <<"  " <<(Pot_tot+Kin_tot)/Num_Of_Atoms ;
		pot_file << "  " <<Sys_Temperature<< "  " << Sys_Pressure << "  " << Xi << "   " << Num_Of_Atoms <<'\n';
	}
}
/***********************************************************************/
/***********************************************************************/
void Allocate_Arrays()
{	
	int i,j;
	Atoms=new Atom_Struct[Max_Num_Of_Atoms];
	for (i=0;i<Max_Num_Of_Atoms;i++) {
		for (j=0;j<3;j++) {
			Atoms[i].R_Vec[j]=0;
			Atoms[i].V_Vec[j]=0;
			Atoms[i].A_Vec[j]=0;
		}
	}
	for(i=0;i<Max_Num_Of_Atoms;i++) {
		(Atoms[i].NeighborList)=new int[max_ng];
		if (Num_Of_Bonds>0)		    Atoms[i].Bonds_Of_Atom=new int[4];
		else						Atoms[i].Bonds_Of_Atom=NULL;
		if (Num_Of_Angles>0)		Atoms[i].Angles_Of_Atom=new int[6];
		else						Atoms[i].Angles_Of_Atom=NULL;
		if (Num_Of_Dihedrals>0)		Atoms[i].Dihedrals_Of_Atom=new int[10];
		else						Atoms[i].Dihedrals_Of_Atom=NULL;
	}

	bondlist=new int [Num_Of_Bonds][3];
	anglelist=new int [Num_Of_Angles][4];
	dihedrallist=new int [Num_Of_Dihedrals][5];
	Delete_Bound_OK=0;
	if ((Boundary_Type.x=='f')|(Boundary_Type.x=='F')) {
		for (i=0;i<Max_Num_Of_Atoms;i++) Atoms[i].Delete_List=0;
	} else if ((Boundary_Type.y=='f')|(Boundary_Type.y=='F')) {
		for (i=0;i<Max_Num_Of_Atoms;i++) Atoms[i].Delete_List=0;
	}else if ((Boundary_Type.z=='f')|(Boundary_Type.z=='F')){
		for (i=0;i<Max_Num_Of_Atoms;i++) Atoms[i].Delete_List=0;
	}
}
/***********************************************************************/
/***********************************************************************/
void Initialization()
{
	double s_T;
	int tp;
	switch (Trj_INIT_MODE) {
		case 'I': {			// Initializing
			srand(time(0));
			for (int i=0;i<Num_Of_Atoms;i++) {
				tp=Atoms[i].tp;
				s_T=sqrt(Boltzmann_const*Target_Temp/AtomMass[tp-1]);
				Atoms[i].V_Vec[0]=Normal_Random(0,s_T);
				Atoms[i].V_Vec[1]=Normal_Random(0,s_T);
				Atoms[i].V_Vec[2]=Normal_Random(0,s_T);
			}
			Calc_Properties_Mid();
			Calc_Properties();
			Immobalize_System();
		} 
            break;
        case 'R': 	Read__Trajectory_Restart(); // read from arestart file.
             break;
        default : cout << "Error Input Mode!" << '\n';
    }
	Fix_Velocity() ;

}
/***********************************************************************/
/***********************************************************************/
void	Get_Input()
{	
	int i,iii;
	double Eps,Sig,Mas,Electric_Charge;
	char D[40], D2;
	InputFile  >> D >> D  >> Trj_INIT_MODE;
	InputFile  >> D >> D  >> Smp_INIT_MODE;
	InputFile >> D  >> D  >>  D ;
	InputFile >> Num_Of_Atoms;
	InputFile >> D  >> D  >>  D ;
	InputFile >> Mem_Increase_factor;
	Max_Num_Of_Atoms=int(floor(Num_Of_Atoms*Mem_Increase_factor));
	InputFile >> D  >> D  >>  D >> D ;
	InputFile >> Num_Of_Atom_Types;
	/***************************************************************************/
	if (StepNum==0) {
		Sigma=new double [Num_Of_Atom_Types];
		Epsilon=new double [Num_Of_Atom_Types];
		AtomMass=new double [Num_Of_Atom_Types];
		q = new double [Num_Of_Atom_Types];
		Sigma_MAT=new double* [Num_Of_Atom_Types];
		for( i=0;i<Num_Of_Atom_Types;i++) *(Sigma_MAT+i)=new double[Num_Of_Atom_Types];
		Epsilon_MAT=new double* [Num_Of_Atom_Types];
		for( i=0;i<Num_Of_Atom_Types;i++) *(Epsilon_MAT+i)=new double[Num_Of_Atom_Types];
	}
	/***************************************************************************/
	InputFile >> D  >> D  >>  D ;
	InputFile >> Num_Of_Bonds;
	InputFile >> D  >> D  >>  D >> D ;
	InputFile >> Num_Of_Bond_Types;
	if (StepNum==0) {
		Harmonic_Bond_Parameter=new Harmonic_B[Num_Of_Bond_Types];
	}
	/***************************************************************************/
	InputFile >> D  >> D  >>  D ;
	InputFile >> Num_Of_Angles;
	InputFile >> D  >> D  >>  D >> D;
	InputFile >> Num_Of_Angle_Types;
	if (StepNum==0) {
		Harmonic_Angle_Parameter=new Harmonic_B[Num_Of_Angle_Types];
	}
	/***************************************************************************/
	InputFile >> D  >> D  >>  D ;
	InputFile >> Num_Of_Dihedrals;
	InputFile >> D  >> D  >>  D >> D;
	InputFile >> Num_Of_Dihedral_Types;
	if (StepNum==0) {
		Harmonic_Dihedral_Parameter=new DIHEDRAL[Num_Of_Dihedral_Types];
	}
	/***************************************************************************/
	InputFile >>	 D  >>	D >>	D >>	D >> Num_Of_Thrust_Sets;
	if (StepNum==0) {
		Thrust_Set = new Propeller[Num_Of_Thrust_Sets];
	}
	/***************************************************************************/
	InputFile >>	 D  >>	D >>	D >>	D >> Num_Of_Vel_Sets;
	if (StepNum==0) {
		Vel_Set=new V_Set[Num_Of_Vel_Sets];
	}
	/***************************************************************************/
	InputFile >>	 D  >>	D >>	D >>	D >> Num_Of_Cons_Temp_Sets;
	if (StepNum==0) {
		cons_Temp_Wall=new C_T_W[Num_Of_Cons_Temp_Sets];
	}
	/***************************************************************************/
	InputFile >>	 D  >>	D >>	D >>	D >> Num_Of_Const_Flux_Sets;
	if (StepNum==0) {
		cons_Flux_Wall=new C_F_W[Num_Of_Const_Flux_Sets]; 
	}
	/***************************************************************************/
	InputFile >>	 D  >>	D >>	D >>	D >> Num_Of_Force_Sets;
	if (StepNum==0) {
		F_Set=new F_Struct[Num_Of_Force_Sets];
	}
	/***************************************************************************/
	InputFile >>	 D  >>	D >>	D >>	D >> Num_Of_Solid_Sets;
	if (StepNum==0) {
		Solid_Wall=new Solid[Num_Of_Solid_Sets];
	}
	/***************************************************************************/
	InputFile >> D  >> D >> D >> D  >> D  >>  D  >> D >> D >> D;
	Coulomb_OK=0;
	for (i=0;i<Num_Of_Atom_Types;i++) {
		InputFile >> iii  >> D >> Eps  >>  Sig >> Mas >> Electric_Charge;
		Epsilon[i]=Eps;
		Sigma[i]=Sig;
		AtomMass[i]=Mas;
		q[i]=Electric_Charge;
		if (fabs(q[i]) > 0.000001)    Coulomb_OK=1;
	}

	sigma_max=Sigma[0];
	sigma_min=Sigma[0];
	for (i=1;i<Num_Of_Atom_Types;i++) {
		if (Sigma[i]>sigma_max) sigma_max=Sigma[i];
		if (Sigma[i]<sigma_min) sigma_min=Sigma[i];
	}
	/****************************************************************************/
	InputFile >> D  >> D  >>  D  >> D >> D >> D >> D;
	for (i=0;i<Num_Of_Bond_Types;i++) {
		InputFile >> iii >> D >> Harmonic_Bond_Parameter[i].K >>Harmonic_Bond_Parameter[i].r0  ;
	}
	InputFile >> D  >> D  >>  D  >> D >> D >> D >> D;
	for (i=0;i<Num_Of_Angle_Types;i++) {
		InputFile >> iii >> D >> Harmonic_Angle_Parameter[i].K>> Harmonic_Angle_Parameter[i].r0  ;
		Harmonic_Angle_Parameter[i].r0*=pi_N/180;
		Harmonic_Angle_Parameter[i].K*=(180*180)/(pi_N*pi_N);

	}
	InputFile >> D  >> D  >>  D  >> D >> D >> D >> D >> D;
	for (i=0; i<Num_Of_Dihedral_Types; i++){
		InputFile >> iii >> D >> Harmonic_Dihedral_Parameter[i].K >> Harmonic_Dihedral_Parameter[i].sign >> Harmonic_Dihedral_Parameter[i].multiplicity;
	}
	InputFile >> D  >> D >>D;
	InputFile >> DeltaT;
	InputFile >>	 D  >>	MaxStep  >> D >> Trj_Write_Step>> D >> Prop_Write_Step>> D >> Restart_Write_Step;
	InputFile >>	 D  >>	D  >> D >> R_Cut_Off_R;
	InputFile >>	 D  >>	D  >> D >> Delta_R_R;
	InputFile >>	 D  >>	R_Cut_Coul;
	R_Cut_Off=R_Cut_Off_R*sigma_max;
	Delta_R=Delta_R_R*sigma_max;
	max_R_Cut=R_Cut_Off+Delta_R;
	if ((R_Cut_Coul>max_R_Cut) && (Coulomb_OK==1)) 	max_R_Cut=R_Cut_Coul;

	InputFile >>	 D  >>	D;
	InputFile >>	 D  >>	X_Low;
	InputFile >>	 D  >>	X_High;
	InputFile >>	 D  >>	Y_Low;
	InputFile >>	 D  >>	Y_High;
	InputFile >>	 D  >>	Z_Low;
	InputFile >>	 D  >>	Z_High;
	InputFile >>	 D  >>	D;
	InputFile >>	 D  >>	X_L_Measure;
	InputFile >>	 D  >>	X_H_Measure;
	InputFile >>	 D  >>	Y_L_Measure;
	InputFile >>	 D  >>	Y_H_Measure;
	InputFile >>	 D  >>	Z_L_Measure;
	InputFile >>	 D  >>	Z_H_Measure;
	InputFile >> D  >> D >> D   >> D >> N_Sam_CELL[0] >> N_Sam_CELL[1]>> N_Sam_CELL[2] ;
	InputFile >>	 D  >> D	>> Sampling_Step >> D  >> Sampling_start;
	InputFile >>	 D  >> D >> Display_Scale;
	InputFile >>	 D  >> D;
	for (i=0;i<3;i++) {
		InputFile >>	 D2 ;
		if ((D2=='x')|(D2=='X')) InputFile >> Boundary_Type.x;
		if ((D2=='y')|(D2=='Y')) InputFile >> Boundary_Type.y;
		if ((D2=='z')|(D2=='Z')) InputFile >> Boundary_Type.z;
	}

	InputFile >> D  >> D >> D >> D >> D >> D;
	for (i = 0; i < Num_Of_Thrust_Sets; i++) {
		InputFile >>	 Thrust_Set[i].Atom_No  >>  Thrust_Set[i].F >> Thrust_Set[i].Atom_No_Direction_1 >> Thrust_Set[i].Atom_No_Direction_2;
	}

	InputFile >> D  >> D >> D >> D >> D >>D;
	for (i=0;i<Num_Of_Vel_Sets;i++) {
		InputFile >>	 Vel_Set[i].tp  >>	Vel_Set[i].vx >>	Vel_Set[i].vy >>	Vel_Set[i].vz ;
	}

	InputFile >> D  >> D >> D >> D >> D >> D >> D>> D;
	for (i=0;i<Num_Of_Cons_Temp_Sets;i++) {
		InputFile >>	 cons_Temp_Wall[i].tp  >>	cons_Temp_Wall[i].vx >>	  cons_Temp_Wall[i].vy >>	  cons_Temp_Wall[i].vz >>   cons_Temp_Wall[i].Temp >>   cons_Temp_Wall[i].K ;
	    cons_Temp_Wall[i].xi=0.0;
	}
	//////////////////////////////////////////
	InputFile >> D  >> D >> D >> D >> D >> D >> D>> D;      
	for (i=0;i<Num_Of_Const_Flux_Sets;i++) {
		InputFile >>	 cons_Flux_Wall[i].tp  >>	cons_Flux_Wall[i].vx >>	  cons_Flux_Wall[i].vy >>	  cons_Flux_Wall[i].vz >>   cons_Flux_Wall[i].Flux >>   cons_Flux_Wall[i].K ;
	}
	//////////////////////////////////////////
	InputFile >> D  >> D >> D >> D >>D;
	for (i=0;i<Num_Of_Force_Sets;i++) {
		InputFile >>	 F_Set[i].tp  >>	F_Set[i].dr >>	F_Set[i].mag;
	}
	InputFile >> D  >> D >> D >> D >> D  >> D >> D  >> D >> D >>D;
	for (i=0;i<Num_Of_Solid_Sets;i++) {
		InputFile  >> Solid_Wall[i].tp >> Solid_Wall[i].init_Vx >> Solid_Wall[i].init_Vy>> Solid_Wall[i].init_Vz >> Solid_Wall[i].fx >> Solid_Wall[i].fy >> Solid_Wall[i].fz >> Solid_Wall[i].si ;
	}
	InputFile >> D  >> D >>  D >> D ;
	InputFile >> Target_Temp >> Thermo_tp>> T_Coupling;
	
	InputFile >> D  >> D >> D >> D >> D ;
	InputFile >> Num_of_Inlets;
	Inlet=new Inlet_Struct[Num_of_Inlets];
	InputFile >> D  >> D >> D >> D >> D >> D >> D >> D >> D  >> D >> D >> D;
	for(i=0;i<Num_of_Inlets;i++) {
		InputFile >> Inlet[i].atom_tp >> Inlet[i].Vx >> Inlet[i].Vy >> Inlet[i].Vz >> Inlet[i].temp>> Inlet[i].Rho;
		InputFile >> Inlet[i].x1 >> Inlet[i].x2 >> Inlet[i].y1 >> Inlet[i].y2 >> Inlet[i].z1 >> Inlet[i].z2;
	}
	InputFile >> D  >> D >> D ;
	InputFile >> Num_Of_Outlets;
	Outlet=new Outlet_Struct[Num_Of_Outlets];
	InputFile >> D  >> D >> D >> D >> D >> D >> D;
	for(i=0;i<Num_Of_Outlets;i++) {
		InputFile >>  Outlet[i].dir >>  Outlet[i].x1 >> Outlet[i].x2 >> Outlet[i].y1 >> Outlet[i].y2 >> Outlet[i].z1 >> Outlet[i].z2;
		if (Outlet[i].dir=='R') {
			Outlet[i].x1=Outlet[i].x2-R_Cut_Off;
		}else if (Outlet[i].dir=='L') {
			Outlet[i].x2=Outlet[i].x1+R_Cut_Off;
		}else if (Outlet[i].dir=='U') {
			Outlet[i].y1=Outlet[i].y2-R_Cut_Off;
		}else if (Outlet[i].dir=='D') {
			Outlet[i].y2=Outlet[i].y1+R_Cut_Off;
		}else if (Outlet[i].dir=='F') {
			Outlet[i].z1=Outlet[i].z2-R_Cut_Off;
		}else if (Outlet[i].dir=='B') {
			Outlet[i].z2=Outlet[i].z1+R_Cut_Off;
		}
	}
		
	X_H_L=X_High-X_Low;
	Y_H_L=Y_High-Y_Low;
	Z_H_L=Z_High-Z_Low;
	ofstream chk("Check.txt");
	chk << " Target_Temp      Thermo_tp    T_Coupling    " << '\n'; 
	chk <<  Target_Temp  <<"                    " <<   Thermo_tp <<"              " <<T_Coupling    ;
	chk.close();
}
/***********************************************************************/
/***********************************************************************/	
void	Get_Data()
{	
	int i,iii;
	int atm1,atm2,atm3, atm4;
	char D[40];
	Allocate_Arrays();
	InputFile2 >> D >>D  >> D  >>  D  >> D >> D ;
	for (i=0;i<Num_Of_Atoms;i++) {
		InputFile2 >> iii >> Atoms[i].R_Vec[0] >> Atoms[i].R_Vec[1] >> Atoms[i].R_Vec[2] >> Atoms[i].tp;
	}
	for (i=0;i<Num_Of_Atoms;i++) {
		Atoms[i].Num_Of_Bonds_Of_Atom=0;
		Atoms[i].Num_Of_Angles_Of_Atom=0;
		Atoms[i].Num_Of_Dihedrals_Of_Atom=0;
	}

	InputFile2 >>  D  >> D  >>  D  >> D  >> D;

	for (i=0;i<Num_Of_Bonds;i++) {
		InputFile2 >> iii >> atm1  >> atm2 >> bondlist[i][2];

		bondlist[i][0]=atm1-1;
		bondlist[i][1]=atm2-1;
		
		Atoms[atm1-1].Bonds_Of_Atom[Atoms[atm1-1].Num_Of_Bonds_Of_Atom]=i;
		Atoms[atm2-1].Bonds_Of_Atom[Atoms[atm2-1].Num_Of_Bonds_Of_Atom]=i;
		Atoms[atm1-1].Num_Of_Bonds_Of_Atom++;
		Atoms[atm2-1].Num_Of_Bonds_Of_Atom++;
	}

	InputFile2 >>  D  >> D  >>  D  >> D  >> D >> D;
	for (i=0;i<Num_Of_Angles;i++) {
		InputFile2 >> iii >> atm1 >> atm2 >> atm3 >> anglelist[i][3];
		anglelist[i][0]=atm1-1;
		anglelist[i][1]=atm2-1;
		anglelist[i][2]=atm3-1;
		Atoms[atm1-1].Angles_Of_Atom[Atoms[atm1-1].Num_Of_Angles_Of_Atom]=i;
		Atoms[atm2-1].Angles_Of_Atom[Atoms[atm2-1].Num_Of_Angles_Of_Atom]=i;
		Atoms[atm3-1].Angles_Of_Atom[Atoms[atm3-1].Num_Of_Angles_Of_Atom]=i;
		Atoms[atm1-1].Num_Of_Angles_Of_Atom++;
		Atoms[atm2-1].Num_Of_Angles_Of_Atom++;
		Atoms[atm3-1].Num_Of_Angles_Of_Atom++;
	}

	InputFile2 >> D  >> D  >>  D  >> D  >> D >> D  >> D;
	for (i=0; i<Num_Of_Dihedrals; i++) {
		InputFile2 >> iii >> atm1 >> atm2 >> atm3 >> atm4 >> dihedrallist[i][4];
		dihedrallist[i][0] = atm1 - 1;
		dihedrallist[i][1] = atm2 - 1;
		dihedrallist[i][2] = atm3 - 1;
		dihedrallist[i][3] = atm4 - 1;
		Atoms[atm1-1].Dihedrals_Of_Atom[Atoms[atm1-1].Num_Of_Dihedrals_Of_Atom] = i;
		Atoms[atm2-1].Dihedrals_Of_Atom[Atoms[atm2-1].Num_Of_Dihedrals_Of_Atom] = i;
		Atoms[atm3-1].Dihedrals_Of_Atom[Atoms[atm3-1].Num_Of_Dihedrals_Of_Atom] = i;
		Atoms[atm4-1].Dihedrals_Of_Atom[Atoms[atm4-1].Num_Of_Dihedrals_Of_Atom] = i;
		Atoms[atm1-1].Num_Of_Dihedrals_Of_Atom++;
		Atoms[atm2-1].Num_Of_Dihedrals_Of_Atom++;
		Atoms[atm3-1].Num_Of_Dihedrals_Of_Atom++;
		Atoms[atm4-1].Num_Of_Dihedrals_Of_Atom++;
	}
	for (i=0;i<Num_Of_Atoms;i++) 	M_system+=AtomMass[Atoms[i].tp-1];
}
/***********************************************************************/
/***********************************************************************/
void 	Set_Parameters()
{
	int i,j;
	for( i=0;i<Num_Of_Atom_Types;i++) {
		for( j=0;j<Num_Of_Atom_Types;j++){
			Sigma_MAT[i][j]=(Sigma[i]+Sigma[j])/2;
			Epsilon_MAT[i][j]=sqrt(Epsilon[i]*Epsilon[j]);
		}
	}
}
/***********************************************************************/
/***********************************************************************/
void X_Periodic_Boundary(int ii) 
{
	if (Atoms[ii].R_Vec[0]<X_Low) {
		Atoms[ii].R_Vec[0]=X_H_L+Atoms[ii].R_Vec[0];
	} else	if (Atoms[ii].R_Vec[0]>X_High) {
		Atoms[ii].R_Vec[0]=-X_H_L+Atoms[ii].R_Vec[0];
	}
}
/***********************************************************************/
/***********************************************************************/
void X_Free_Boundary(int ii) 
{
	double max_R=sigma_max*R_Cut_Off_R;
	if ((Atoms[ii].R_Vec[0]<X_Low_E) | (Atoms[ii].R_Vec[0]>X_High_E)) {
		Atoms[ii].Delete_List=1;
		Delete_Bound_OK=1;
	}
}
/***********************************************************************/
/***********************************************************************/
void Y_Periodic_Boundary(int ii) 
{
	if (Atoms[ii].R_Vec[1]<Y_Low) {
		Atoms[ii].R_Vec[1]=Y_H_L+Atoms[ii].R_Vec[1];
	}else 	if (Atoms[ii].R_Vec[1]>Y_High) {
		Atoms[ii].R_Vec[1]=-Y_H_L+Atoms[ii].R_Vec[1];
	}
}
/***********************************************************************/
/***********************************************************************/
void Y_Free_Boundary(int ii) 
{
	double max_R=sigma_max*R_Cut_Off_R;
	if ((Atoms[ii].R_Vec[1]<Y_Low_E) | (Atoms[ii].R_Vec[1]>Y_High_E)) {
		Atoms[ii].Delete_List=1;
		Delete_Bound_OK=1;
	}
}
/***********************************************************************/
/***********************************************************************/
void Z_Periodic_Boundary(int ii) 
{
	if (Atoms[ii].R_Vec[2]<Z_Low) {
		Atoms[ii].R_Vec[2]=Z_H_L+Atoms[ii].R_Vec[2];
	}else 	if (Atoms[ii].R_Vec[2]>Z_High) {
		Atoms[ii].R_Vec[2]=-Z_H_L+Atoms[ii].R_Vec[2];
	}
}
/***********************************************************************/
/***********************************************************************/
void Z_Free_Boundary(int ii) 
{
	double max_R=sigma_max*R_Cut_Off_R;
	if ((Atoms[ii].R_Vec[2]<Z_Low_E) | (Atoms[ii].R_Vec[2]>Z_High_E)) {
		Atoms[ii].Delete_List=1;
		Delete_Bound_OK=1;
	}
}
/***********************************************************************/
/***********************************************************************/
void Dump_Trajectory() 
{
	int 	STP_I;
	char STP_C[10],STR1[20]="Trajectory/TRJ.";
	double XL,XH,YL,YH,ZL,ZH;
	int i;
	STP_I = StepNum;
	itoa(STP_I,STP_C,10);
	strcat(STR1, STP_C);
	ofstream TJC(STR1);
	
	XL=X_Low-0.1*(X_High-X_Low);
	XH=X_High+0.1*(X_High-X_Low);
	YL=Y_Low-0.1*(Y_High-Y_Low);
	YH=Y_High+0.1*(Y_High-Y_Low);
	ZL=Z_Low-0.1*(Z_High-Z_Low);
	ZH=Z_High+0.1*(Z_High-Z_Low);
	cout << "Step=" << StepNum << "          Num="<<  Num_Of_Atoms << '\n';
	TJC << "ITEM: TIMESTEP" << '\n';
	TJC << StepNum << '\n';
	TJC << "ITEM: NUMBER OF ATOMS" << '\n';
	TJC << Num_Of_Atoms<< '\n';
	TJC << "ITEM: BOX BOUNDS" << '\n';
	TJC << XL << "  "  << XH << '\n';
	TJC << YL << "  "  << YH << '\n';
	TJC << ZL << "  "  << ZH << '\n';
	TJC << "ITEM: ATOMS" << '\n';
	for(i=0;i<Num_Of_Atoms;i++) {
		TJC << i+1 << "  " << Atoms[i].tp<< "  " << Atoms[i].R_Vec[0]/Display_Scale << "  " << Atoms[i].R_Vec[1]/Display_Scale << "  " << Atoms[i].R_Vec[2]/Display_Scale << '\n' ;
	}
	TJC.close();
}
/***********************************************************************/
/***********************************************************************/
void External_Force()
{
	int i;
	for (i=0;i<Num_Of_Force_Sets;i++) {
		Add_Force(F_Set[i].tp,F_Set[i].dr,F_Set[i].mag);
	}
}
/***********************************************************************/
/***********************************************************************/
void Add_Force(int Atom_t,char direction,double Ext_Force) 
{
	int i,j;
	if ((direction=='x')|(direction=='X')) j=0;
	if ((direction=='y')|(direction=='Y')) j=1;
	if ((direction=='z')|(direction=='Z')) j=2;
	for (i=0;i<Num_Of_Atoms;i++){
		if (Atoms[i].tp==Atom_t) Atoms[i].force[j]+=Ext_Force;
	}
}		
/***********************************************************************/
/***********************************************************************/
void Create_Neighbourhood() 
{ 
	int i,i1,i2,i3,i4,i5,i6,j,k3,k4,Inc,i7;
	int nmn=(2*MMM+1)*(2*MMM+1)*(2*MMM+1);
	double x_ij,y_ij,z_ij,r_ij;
	double sigma;
	New_Cells();
	for (i6=0;i6<Num_Of_Atoms;i6++) Atoms[i6].NeighborList[0]=0;
	for (i5=0;i5<N_CELLS;i5++) {
		for (k3=0;k3<nmn;k3++) {
			k4=Neighb_Cell[i5][k3];
			Considered_Cell[k4]=0;
		}
		for (k3=0;k3<nmn;k3++) {
			k4=Neighb_Cell[i5][k3];
			if (Considered_Cell[k4]==0) {
				Considered_Cell[k4]=1;
				for (i6=1;i6<CELL_MEMBER[i5][0]+1;i6++) {
					i=CELL_MEMBER[i5][i6];
					for (i2=1;i2<CELL_MEMBER[k4][0]+1;i2++) {
						j=CELL_MEMBER[k4][i2];
						if (j>i) {
						/************************************************/
							//  Excluding pairs of the same bond or angle   
							Inc=1;	
							for (i3=0;i3<Atoms[i].Num_Of_Bonds_Of_Atom;i3++) {
								for (i4=0;i4<Atoms[j].Num_Of_Bonds_Of_Atom;i4++) {
									if (Atoms[i].Bonds_Of_Atom[i3]==Atoms[j].Bonds_Of_Atom[i4]) Inc=0;
								}
							}
							for (i3=0;i3<Atoms[i].Num_Of_Angles_Of_Atom;i3++) {
								for (i4=0;i4<Atoms[j].Num_Of_Angles_Of_Atom;i4++) {
									if (Atoms[i].Angles_Of_Atom[i3]==Atoms[j].Angles_Of_Atom[i4]) Inc=0;
								}
							}
							for (i3=0;i3<Atoms[i].Num_Of_Dihedrals_Of_Atom;i3++) {
								for (i4=0;i4<Atoms[j].Num_Of_Dihedrals_Of_Atom;i4++) {
									if (Atoms[i].Dihedrals_Of_Atom[i3]==Atoms[j].Dihedrals_Of_Atom[i4]) Inc=0;
								}
							}
							for (i1=0;i1<Num_Of_Vel_Sets;i1++) {
								for (i7=0;i7<Num_Of_Vel_Sets;i7++) {
									if ((Atoms[i].tp==Vel_Set[i1].tp)&&(Atoms[j].tp==Vel_Set[i7].tp)) Inc=0;   // This Inc is necessary because several checks are done here.
								}
							}
/*
							for (i1=0;i1<Num_Of_Cons_Temp_Sets;i1++) {
								for (i7=0;i7<Num_Of_Cons_Temp_Sets;i7++) {
                                    if ((Atoms[i].tp==cons_Temp_Wall[i1].tp)&&(Atoms[j].tp==cons_Temp_Wall[i7].tp)) Inc=1;
								}
							}
*/
							for (i7=0;i7<Num_Of_Solid_Sets;i7++) {
								if ((Atoms[i].tp==Solid_Wall[i7].tp)&&(Atoms[j].tp==Solid_Wall[i7].tp)) Inc=0;   // This Inc is necessary because several checks are done here.
							}
							if (Inc==1) {
								z_ij=(Atoms[i].R_Vec[2]-Atoms[j].R_Vec[2]);
								y_ij=(Atoms[i].R_Vec[1]-Atoms[j].R_Vec[1]);
								x_ij=(Atoms[i].R_Vec[0]-Atoms[j].R_Vec[0]);
								if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
									if ((x_ij)>(X_H_L)/2) x_ij=x_ij-X_H_L;
									else if ((x_ij)<-X_H_L/2) x_ij=x_ij+X_H_L;
								}
								if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
									if ((y_ij)>(Y_H_L)/2) y_ij=y_ij-Y_H_L;
									else if ((y_ij)<-Y_H_L/2) y_ij=y_ij+Y_H_L;
								}
								
								if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
									if ((z_ij)>Z_H_L/2) z_ij=z_ij-Z_H_L;
									else if ((z_ij)<-Z_H_L/2) z_ij=z_ij+Z_H_L;
								}	
								
								r_ij=sqrt(x_ij*x_ij+y_ij*y_ij+z_ij*z_ij);
								sigma=Sigma_MAT[Atoms[i].tp-1][Atoms[j].tp-1];
								R_Cut_Off=R_Cut_Off_R*sigma;
								Delta_R=Delta_R_R*sigma;
								max_R_Cut=R_Cut_Off+Delta_R;
								if ((R_Cut_Coul>max_R_Cut) && (fabs(q[Atoms[i].tp-1])>0.000001) ) 	max_R_Cut=R_Cut_Coul;					
								if ((r_ij<max_R_Cut)) { 
									Atoms[i].NeighborList[0]++;
									Atoms[i].NeighborList[Atoms[i].NeighborList[0]]=j;
								}
							}
						}
					}
				}
			} 
		}
	}
}
/***********************************************************************/
/***********************************************************************/
void Create_Neighbourhood_of_Atom(int i) 
{ 
	int i1,i2,i3,i4,i5,j,k3,k4,Inc,i7;
	int nmn=(2*MMM+1)*(2*MMM+1)*(2*MMM+1);
	double x_ij,y_ij,z_ij,r_ij;
	double sigma;
	Cell_Det(i);
	i5=Atoms[i].Cell;
	Atoms[i].NeighborList[0]=0;
	for (k3=0;k3<nmn;k3++) {
		k4=Neighb_Cell[i5][k3];
		Considered_Cell[k4]=0;
	
	}
	for (k3=0;k3<nmn;k3++) {
		k4=Neighb_Cell[i5][k3];
		if (Considered_Cell[k4]==0) {
			Considered_Cell[k4]=1;
		
			for (i2=1;i2<CELL_MEMBER[k4][0]+1;i2++) {
				j=CELL_MEMBER[k4][i2];
				if (j!=i) {
				/************************************************/
					//  Excluding pairs of the same bond or angle   
					Inc=1;	
					for (i3=0;i3<Atoms[i].Num_Of_Bonds_Of_Atom;i3++) {
						for (i4=0;i4<Atoms[j].Num_Of_Bonds_Of_Atom;i4++) {
							if (Atoms[i].Bonds_Of_Atom[i3]==Atoms[j].Bonds_Of_Atom[i4]) Inc=0;
						}
					}
					for (i3=0;i3<Atoms[i].Num_Of_Angles_Of_Atom;i3++) {
						for (i4=0;i4<Atoms[j].Num_Of_Angles_Of_Atom;i4++) {
							if (Atoms[i].Angles_Of_Atom[i3]==Atoms[j].Angles_Of_Atom[i4]) Inc=0;
						}
					}
					for (i1=0;i1<Num_Of_Vel_Sets;i1++) {
						for (i7=0;i7<Num_Of_Vel_Sets;i7++) {
							if ((Atoms[i].tp==Vel_Set[i1].tp)&&(Atoms[j].tp==Vel_Set[i7].tp)) Inc=0;   // This Inc is necessary because several checks are done here.
						}
					}
/*		               
					for (i1=0;i1<Num_Of_Cons_Temp_Sets;i1++) {
						for (i7=0;i7<Num_Of_Cons_Temp_Sets;i7++) {
                                  if ((Atoms[i].tp==cons_Temp_Wall[i1].tp)&&(Atoms[j].tp==cons_Temp_Wall[i7].tp)) Inc=1;
						}
					}
*/								
					for (i7=0;i7<Num_Of_Solid_Sets;i7++) {
						if ((Atoms[i].tp==Solid_Wall[i7].tp)&&(Atoms[j].tp==Solid_Wall[i7].tp)) Inc=0;   // This Inc is necessary because several checks are done here.
					}
					if (Inc==1) {
						z_ij=(Atoms[i].R_Vec[2]-Atoms[j].R_Vec[2]);
						y_ij=(Atoms[i].R_Vec[1]-Atoms[j].R_Vec[1]);
						x_ij=(Atoms[i].R_Vec[0]-Atoms[j].R_Vec[0]);
						if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
							if ((x_ij)>(X_H_L)/2) x_ij=x_ij-X_H_L;
							else if ((x_ij)<-X_H_L/2) x_ij=x_ij+X_H_L;
						}
						if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
							if ((y_ij)>(Y_H_L)/2) y_ij=y_ij-Y_H_L;
							else if ((y_ij)<-Y_H_L/2) y_ij=y_ij+Y_H_L;
						}
						
						if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
							if ((z_ij)>Z_H_L/2) z_ij=z_ij-Z_H_L;
							else if ((z_ij)<-Z_H_L/2) z_ij=z_ij+Z_H_L;
						}	
						
						r_ij=sqrt(x_ij*x_ij+y_ij*y_ij+z_ij*z_ij);
						sigma=Sigma_MAT[Atoms[i].tp-1][Atoms[j].tp-1];
						R_Cut_Off=R_Cut_Off_R*sigma;
						Delta_R=Delta_R_R*sigma;
						max_R_Cut=R_Cut_Off+Delta_R;
						if ((R_Cut_Coul>max_R_Cut) && (fabs(q[Atoms[i].tp-1])>0.000001) ) 	max_R_Cut=R_Cut_Coul;					
						if ((r_ij<max_R_Cut)) { 
							Atoms[i].NeighborList[0]++;
							Atoms[i].NeighborList[Atoms[i].NeighborList[0]]=j;
						}
					}
				}
			}
	
		} 
	}
	
}
/***********************************************************************/
/***********************************************************************/
void Make_Solid_List () 
{
	int i,kk[1000],j,j1,j2;
	double Xc,Yc,Zc,Rx,Ry,Rz;
	for (i=0;i<1000;i++) kk[i]=-1;
	for (j=0;j<Num_Of_Solid_Sets;j++) {
		Solid_Wall[j].Num=0;
	}
	for (i=0;i<Num_Of_Atoms;i++) {
		for (j=0;j<Num_Of_Solid_Sets;j++) {
			if (Atoms[i].tp==Solid_Wall[j].tp) Solid_Wall[j].Num++;
			if (Atoms[i].tp==Solid_Wall[j].tp)  break;
		}
	}
	Solid_Atoms=new int* [Num_Of_Solid_Sets];
	for( i=0;i<Num_Of_Solid_Sets;i++) *(Solid_Atoms+i)=new int[Solid_Wall[i].Num]; 

	for (i=0;i<Num_Of_Atoms;i++) {
		for (j=0;j<Num_Of_Solid_Sets;j++) {
			if (Atoms[i].tp==Solid_Wall[j].tp) {
				kk[j]++;
				Solid_Atoms[j][kk[j]]=i;
			}
			if (Atoms[i].tp==Solid_Wall[j].tp)  break;
		}
	}
	for (j=0;j<Num_Of_Solid_Sets;j++) {
		for (j1=0;j1<Solid_Wall[j].Num;j1++) {
			j2=Solid_Atoms[j][j1];
			Atoms[j2].V_Vec[0]=Solid_Wall[j].init_Vx;
			Atoms[j2].V_Vec[1]=Solid_Wall[j].init_Vy;
			Atoms[j2].V_Vec[2]=Solid_Wall[j].init_Vz;
		}
	}

	for (j=0;j<Num_Of_Solid_Sets;j++) {
		Solid_Wall[j].I=0;
		Xc=0;
		Yc=0;
		Zc=0;
		for (j1=0;j1<Solid_Wall[j].Num;j1++) {
			j2=Solid_Atoms[j][j1];
			Xc+=Atoms[j2].R_Vec[0];
			Yc+=Atoms[j2].R_Vec[1];
			Zc+=Atoms[j2].R_Vec[2];
		}
		Xc=Xc/Solid_Wall[j].Num;
		Yc=Yc/Solid_Wall[j].Num;
		Zc=Zc/Solid_Wall[j].Num;
		for (j1=0;j1<Solid_Wall[j].Num;j1++) {
			j2=Solid_Atoms[j][j1];
			Rx=Atoms[j2].R_Vec[0]-Xc;
			Ry=Atoms[j2].R_Vec[1]-Yc;
			Rz=Atoms[j2].R_Vec[2]-Zc;

			Solid_Wall[j].I+=Rx*Rx+Ry*Ry+Rz*Rz;			
		}
		Solid_Wall[j].I=Solid_Wall[j].I*AtomMass[Solid_Wall[j].tp-1];
	}
}
/***********************************************************************/
/***********************************************************************/
void Fix_Velocity() 
{
	int i,j;
	for (i=0;i<N_freezed_mols;i++) {
		for (j=0;j<Num_Of_Vel_Sets;j++) {
			if (Atoms[freezed_mols[i]].tp==Vel_Set[j].tp) {
				Atoms[freezed_mols[i]].V_Vec[0]=Vel_Set[j].vx;
				Atoms[freezed_mols[i]].V_Vec[1]=Vel_Set[j].vy;
				Atoms[freezed_mols[i]].V_Vec[2]=Vel_Set[j].vz;
				Atoms[freezed_mols[i]].A_Vec[0]=0;
				Atoms[freezed_mols[i]].A_Vec[1]=0;
				Atoms[freezed_mols[i]].A_Vec[2]=0;
			}
			if (Atoms[freezed_mols[i]].tp==Vel_Set[j].tp) break;
		}
	}
}
/***********************************************************************/
/***********************************************************************/
void Wall_Temp_Scaling() 
{
	int i,j,I;
	double U,V,W;
	for (i=0;i<Num_Of_Cons_Temp_Sets;i++) {
		U=cons_Temp_Wall[i].vx;
		V=cons_Temp_Wall[i].vy;
		W=cons_Temp_Wall[i].vz;
		Kin_tot_Measure_Wall[i]=0;
	    for (j=0;j<(Cons_Temp_Atoms[i][0]);j++) {
			I=Cons_Temp_Atoms[i][j+1];
	        Kin_tot_Measure_Wall[i]+=(Atoms[I].V_Vec[0]-U)*(Atoms[I].V_Vec[0]-U)+(Atoms[I].V_Vec[1]-V)*(Atoms[I].V_Vec[1]-V)+(Atoms[I].V_Vec[2]-W)*(Atoms[I].V_Vec[2]-W);//+(Sim_Field[I][5]-W)*(Sim_Field[I][5]-W);
		}
		Wall_Temp[i]=(1./3)*Kin_tot_Measure_Wall[i]*(AtomMass[cons_Temp_Wall[i].tp-1]/(Boltzmann_const*Cons_Temp_Atoms[i][0]));
	}
	for (i=0;i<Num_Of_Cons_Temp_Sets;i++) {
		U=cons_Temp_Wall[i].vx;
		V=cons_Temp_Wall[i].vy;
		W=cons_Temp_Wall[i].vz;
	    for (j=0;j<(Cons_Temp_Atoms[i][0]);j++) {
			I=Cons_Temp_Atoms[i][j+1];
            Atoms[I].V_Vec[0]=U+(Atoms[I].V_Vec[0]-U)*sqrt(cons_Temp_Wall[i].Temp/Wall_Temp[i]);
			Atoms[I].V_Vec[1]=V+(Atoms[I].V_Vec[1]-V)*sqrt(cons_Temp_Wall[i].Temp/Wall_Temp[i]);
			Atoms[I].V_Vec[2]=W+(Atoms[I].V_Vec[2]-W)*sqrt(cons_Temp_Wall[i].Temp/Wall_Temp[i]);
		}
	}	    
}
/***********************************************************************/
/***********************************************************************/
void Make_Const_Vel_List () 
{
	int i,kk,j;
	kk=-1;
	N_freezed_mols=0;
	for (i=0;i<Num_Of_Atoms;i++) {
		for (j=0;j<Num_Of_Vel_Sets;j++) {
			if (Atoms[i].tp==Vel_Set[j].tp) N_freezed_mols++;
			if (Atoms[i].tp==Vel_Set[j].tp)  break;
		}
	}
	freezed_mols=new int [N_freezed_mols];
	for (i=0;i<Num_Of_Atoms;i++) {
		for (j=0;j<Num_Of_Vel_Sets;j++) {
			if (Atoms[i].tp==Vel_Set[j].tp) {
				kk++;
				freezed_mols[kk]=i;
			}
			if (Atoms[i].tp==Vel_Set[j].tp)  break;
		}
	}
}
/***********************************************************************/
/***********************************************************************/
void Make_Const_Temp_and_Flux_Wall_List() 
{

	int i,j,J,*kk=NULL,*counter=NULL;
	/************ CONSTANT TEMPERATURE *************/
	kk=new int [Num_Of_Cons_Temp_Sets];
	counter=new int [Num_Of_Cons_Temp_Sets];
	for (j=0;j<(Num_Of_Cons_Temp_Sets);j++) {
		counter[j]=0;
		kk[j]=0;
	}
	for (i=0;i<Num_Of_Atoms;i++) {
		for (j=0;j<(Num_Of_Cons_Temp_Sets);j++) {
			if (Atoms[i].tp==cons_Temp_Wall[j].tp) {
				counter[j]++;
			}
			if (Atoms[i].tp==cons_Temp_Wall[j].tp)  break;
		}
	}
	max_Wall_mol=0;
	for (j=0;j<Num_Of_Cons_Temp_Sets;j++) {
	   if (counter[j]>max_Wall_mol) {
	      max_Wall_mol=counter[j];
	   }
	}
	X_e_T_Wall=new double *[Num_Of_Cons_Temp_Sets];
    Y_e_T_Wall=new double *[Num_Of_Cons_Temp_Sets];
	Z_e_T_Wall=new double *[Num_Of_Cons_Temp_Sets];

	for(i=0;i<Num_Of_Cons_Temp_Sets;i++){
		*(X_e_T_Wall+i)=new double[max_Wall_mol];
		*(Y_e_T_Wall+i)=new double[max_Wall_mol];
		*(Z_e_T_Wall+i)=new double[max_Wall_mol];

	}

	Cons_Temp_Atoms=new int *[Num_Of_Cons_Temp_Sets];
	for(i=0;i<Num_Of_Cons_Temp_Sets;i++){
		*(Cons_Temp_Atoms+i)=new int[max_Wall_mol+1];
	}

	for (i=0;i<Num_Of_Atoms;i++) {
		for (j=0;j<Num_Of_Cons_Temp_Sets;j++) {
			if (Atoms[i].tp==cons_Temp_Wall[j].tp) {
				kk[j]++;
				Cons_Temp_Atoms[j][kk[j]]=i;
                Cons_Temp_Atoms[j][0]=kk[j];
			}
			if (Atoms[i].tp==cons_Temp_Wall[j].tp)  break;
		}
	}
	    for (i=0; i< Num_Of_Cons_Temp_Sets;i++) {
		    for (j=0; j<Cons_Temp_Atoms[i][0];j++) {
				J=Cons_Temp_Atoms[i][j+1];
	           X_e_T_Wall[i][j]=Atoms[J].R_Vec[0];
               Y_e_T_Wall[i][j]=Atoms[J].R_Vec[1];
			   Z_e_T_Wall[i][j]=Atoms[J].R_Vec[2];
		    }
	    }
	Wall_Temp=new double [Num_Of_Cons_Temp_Sets];
	Kin_tot_Measure_Wall=new double [Num_Of_Cons_Temp_Sets];
	/************ CONSTANT FLUX *************/
	int *kk2=NULL,*counter2=NULL;
	kk2=new int [Num_Of_Const_Flux_Sets];
	counter2=new int [Num_Of_Const_Flux_Sets];
	for (j=0;j<(Num_Of_Const_Flux_Sets);j++) {
		counter2[j]=0;
		kk2[j]=0;
	}
	for (i=0;i<Num_Of_Atoms;i++) {
		for (j=0;j<(Num_Of_Const_Flux_Sets);j++) {
			if (Atoms[i].tp==cons_Flux_Wall[j].tp) {
				counter2[j]++;
			}
			if (Atoms[i].tp==cons_Flux_Wall[j].tp)  break;
		}
	}
	max_Wall_mol=0;
	for (j=0;j<Num_Of_Const_Flux_Sets;j++) {
	   if (counter2[j]>max_Wall_mol) {
	      max_Wall_mol=counter2[j];
	   }
	}
	X_e_F_Wall=new double *[Num_Of_Const_Flux_Sets];
    Y_e_F_Wall=new double *[Num_Of_Const_Flux_Sets];
	Z_e_F_Wall=new double *[Num_Of_Const_Flux_Sets];

	for(i=0;i<Num_Of_Const_Flux_Sets;i++){
		*(X_e_F_Wall+i)=new double[max_Wall_mol];
		*(Y_e_F_Wall+i)=new double[max_Wall_mol];
		*(Z_e_F_Wall+i)=new double[max_Wall_mol];
	}

	Const_Flux_Atoms=new int *[Num_Of_Const_Flux_Sets];
	for(i=0;i<Num_Of_Const_Flux_Sets;i++){
		*(Const_Flux_Atoms+i)=new int[max_Wall_mol+1];
	}

	for (i=0;i<Num_Of_Atoms;i++) {
		for (j=0;j<Num_Of_Const_Flux_Sets;j++) {
			if (Atoms[i].tp==cons_Flux_Wall[j].tp) {
				kk2[j]++;
				Const_Flux_Atoms[j][kk2[j]]=i;
                Const_Flux_Atoms[j][0]=kk2[j];
			}
			if (Atoms[i].tp==cons_Flux_Wall[j].tp)  break;
		}
	}
	    for (i=0; i< Num_Of_Const_Flux_Sets;i++) {
		    for (j=0; j<Const_Flux_Atoms[i][0];j++) {
	           X_e_F_Wall[i][j]=Atoms[J].R_Vec[0];
               Y_e_F_Wall[i][j]=Atoms[J].R_Vec[1];
			   Z_e_F_Wall[i][j]=Atoms[J].R_Vec[2];
		    }
	    }
	Wall_Flux=new double [Num_Of_Const_Flux_Sets];
}
/***********************************************************************/
/***********************************************************************/
double Uniform_Random(double x1,double x2)
{ 
	return x1+(x2-x1)*(double(rand()%1000))/999;
}
/***********************************************************************/
/***********************************************************************/
void Make_Cell_Member()
{
	int i1,i2;
	for (i1=0;i1<N_CELLS;i1++) {
		CELL_MEMBER[i1][0]=0;
		for (i2=1;i2<Max_N_Of_Cell;i2++) CELL_MEMBER[i1][i2]=-1;
	}
	for (i1=0;i1<Num_Of_Atoms;i1++) {
		Cell_Det(i1) ;
		CELL_MEMBER[Atoms[i1].Cell][0]++;
		CELL_MEMBER[Atoms[i1].Cell][CELL_MEMBER[Atoms[i1].Cell][0]]=i1;
	}
}
/***********************************************************************/
/***********************************************************************/
void MAKE_Neighb_Cell() 
{
	int KK1,KK2,KK3,i3,i,i1,i2,N1,N2,N3,Nk;
	/******************************************************/
	double L;
	N_CELL[0]=N_CELL[1]=N_CELL[2]=1;
	N_CELLS=1;
	R_Cut_Off=R_Cut_Off_R*sigma_max;
	X_High_E=X_High;   // E : Extended
	X_Low_E=X_Low;
	Y_High_E=Y_High;
	Y_Low_E=Y_Low;
	Z_High_E=Z_High;
	Z_Low_E=Z_Low;
		/*

	if ((Boundary_Type.x=='F') | (Boundary_Type.x=='f') ) 	{
		X_High_E=X_High+R_Cut_Off;
		X_Low_E=X_Low-R_Cut_Off;
	}

	if ((Boundary_Type.y=='f') | (Boundary_Type.y=='F') ) 	{
		Y_High_E=Y_High+R_Cut_Off;
		Y_Low_E=Y_Low-R_Cut_Off;
	}

	if ((Boundary_Type.z=='f') | (Boundary_Type.z=='F') ) 	{
		Z_High_E=Z_High+R_Cut_Off;
		Z_Low_E=Z_Low-R_Cut_Off;
	}
*/

	for (i1=0;i1<Dimension;i1++) {
		if (i1==0) L=X_High_E-X_Low_E;
		else if(i1==1) L=Y_High_E-Y_Low_E;
		else	L=Z_High_E-Z_Low_E;
		R_Cut_Off=R_Cut_Off_R*sigma_max;
		N_CELL[i1]=int(floor(MMM*L/(max_R_Cut)));
		if (N_CELL[i1]==0) N_CELL[i1]=1;		// corrected1
		Dx[i1]=L/N_CELL[i1];
		N_CELLS*=N_CELL[i1];
	}
	Max_N_Of_Cell=int(floor(Rho_max*(Dx[0]*Dx[1]*Dx[2])/(sigma_min*sigma_min*sigma_min)));
	Considered_Cell=new int[N_CELLS];
	Neighb_Cell=new int* [N_CELLS];
	int nmn=(2*MMM+1)*(2*MMM+1)*(2*MMM+1);
	for( i1=0;i1<N_CELLS;i1++) *(Neighb_Cell+i1)=new int[nmn]; 
	CELL_MEMBER=new int* [N_CELLS];
	for( i1=0;i1<N_CELLS;i1++) *(CELL_MEMBER+i1)=new int[Max_N_Of_Cell]; 
	for (i1=0;i1<N_CELLS;i1++) {
		CELL_MEMBER[i1][0]=0;
		for (i2=1;i2<Max_N_Of_Cell;i2++) CELL_MEMBER[i1][i2]=-1;
	}

	for (i1=0;i1<Num_of_Inlets;i1++)   Inlet[i1].N_C=0;
	for (i1=0;i1<Num_Of_Outlets;i1++)   Outlet[i1].N_C=0;

	for (i=0;i<N_CELLS;i++) {
		KK1=i%N_CELL[0];
		KK2=(i%(N_CELL[0]*N_CELL[1]))/N_CELL[0];
		KK3=i/(N_CELL[0]*N_CELL[1]);
		Nk=0;
		for (i1=-MMM;i1<MMM+1;i1++) {
			for (i2=-MMM;i2<MMM+1;i2++) {
				for (i3=-MMM;i3<MMM+1;i3++) {
					N1=(KK1+i1);
					while (N1<0) N1+=N_CELL[0];
					while (N1>N_CELL[0]-1)	N1-=N_CELL[0];

					N2=(KK2+i2);
					while (N2<0) N2+=N_CELL[1];
					while (N2>N_CELL[1]-1)	N2-=N_CELL[1];
	
					N3=(KK3+i3);
					while (N3<0) N3+=N_CELL[2];
					while (N3>N_CELL[2]-1)	N3-=N_CELL[2];

				
					Neighb_Cell[i][Nk]=N_CELL[1]*N_CELL[0]*N3+N_CELL[0]*N2+N1;
					Nk++;
				}
			}
		}

		for (i1=0;i1<Num_of_Inlets;i1++) {
			if (Belong_Cell_to_Inlet(i,i1)==1) {
				Inlet[i1].N_C++;
				Inlet[i1].Cells[Inlet[i1].N_C-1]=i;
				break;
			}
		}
		for (i1=0;i1<Num_Of_Outlets;i1++) {
			if (Belong_Cell_to_Outlet(i,i1)==1) {
				Outlet[i1].N_C++;
				Outlet[i1].Cells[Outlet[i1].N_C-1]=i;
				break;
			}
		}	
	}

	for (i1=0;i1<Num_of_Inlets;i1++) {
		for (i2=0;i2<Inlet[i1].N_C;i2++) {
			Inlet[i1].Cell_N_Res[i2]=0;
		}
	}
}
/***********************************************************************/
/***********************************************************************/
void  Cell_Det(int i1) 
{
int 	KK=0,KK1,KK2,KK3;
		KK1=int(floor((Atoms[i1].R_Vec[0]-X_Low_E)/Dx[0]));
		if (KK1==-1) KK1=0;
		else if (KK1==N_CELL[0]) KK1=N_CELL[0]-1;
		KK2=int(floor((Atoms[i1].R_Vec[1]-Y_Low_E)/Dx[1]));
		if (KK2==-1) KK2=0;
		else if (KK2==N_CELL[1]) KK2=N_CELL[1]-1;
		KK3=int(floor((Atoms[i1].R_Vec[2]-Z_Low_E)/Dx[2]));
		if (KK3==-1) KK3=0;
		else if (KK3==N_CELL[2]) KK3=N_CELL[2]-1;
		Atoms[i1].Cell=N_CELL[0]*N_CELL[1]*KK3+N_CELL[0]*KK2+KK1;
}
/***********************************************************************/
/***********************************************************************/
void Refresh_Cell(int i1,int i2,int i3)
{
	int i;
	for (i=1;i<CELL_MEMBER[i1][0];i++) {
		if (CELL_MEMBER[i1][i]==i3) {
			CELL_MEMBER[i1][i]=CELL_MEMBER[i1][CELL_MEMBER[i1][0]];
			CELL_MEMBER[i1][CELL_MEMBER[i1][0]]=-1;

			break;
		}
	}
	CELL_MEMBER[i1][0]+=-1;
	CELL_MEMBER[i2][0]+=1;
	CELL_MEMBER[i2][CELL_MEMBER[i2][0]]=i3;
}

/***********************************************************************/
/***********************************************************************/
void New_Cells () 
{
	int I,I2,I3;
	for (I=0;I<Num_Of_Atoms;I++) {
		I3=Atoms[I].Cell;
		Cell_Det(I);
		I2=Atoms[I].Cell;
		if (I2!=I3) 	Refresh_Cell(I3,I2,I);
	}
}
/***********************************************************************/
/***********************************************************************/
void Devide_Geo() 
{
	int i,KKK1,KKK2,KKK3;
	
	N_Sam_CELLS=N_Sam_CELL[0]*N_Sam_CELL[1]*N_Sam_CELL[2];
	Dx_s[0]=(X_High_E-X_Low_E)/N_Sam_CELL[0];
	Dx_s[1]=(Y_High_E-Y_Low_E)/N_Sam_CELL[1];
	Dx_s[2]=(Z_High_E-Z_Low_E)/N_Sam_CELL[2];
	Volume_Of_Sam_Cell=Dx_s[0]*Dx_s[1]*Dx_s[2];

	Max_N_Of_Sam_Cell=int(floor(MMM*MMM*MMM*(Max_N_Of_Cell*Volume_Of_Sam_Cell/(max_R_Cut*max_R_Cut*max_R_Cut))))+1;
	Sam_CELL_MEMBER=new int* [N_Sam_CELLS];
	for(int i1=0;i1<N_Sam_CELLS;i1++) *(Sam_CELL_MEMBER+i1)=new int[Max_N_Of_Sam_Cell];

	Vel_Of_Sam_Cell=new double[N_Sam_CELLS][3];
	X_Of_Sam_Cell=new double[N_Sam_CELLS];
	Y_Of_Sam_Cell=new double[N_Sam_CELLS];
	Z_Of_Sam_Cell=new double[N_Sam_CELLS];
	KE_Of_Sam_Cell=new double[N_Sam_CELLS];
	Density_Of_Sam_Cell=new double[N_Sam_CELLS];

	Flux_KE_Of_X_Of_Sam_Cell=new double [N_Sam_CELLS];
	Flux_KE_Of_Y_Of_Sam_Cell=new double [N_Sam_CELLS];
	Flux_KE_Of_Z_Of_Sam_Cell=new double [N_Sam_CELLS];
	Flux_pot_of_Cell_x=new double [N_Sam_CELLS];
	Flux_pot_of_Cell_y=new double [N_Sam_CELLS];
	Flux_pot_of_Cell_z=new double [N_Sam_CELLS];
	Flux_pres_work_Cell_x=new double [N_Sam_CELLS];
	Flux_pres_work_Cell_y=new double [N_Sam_CELLS];
	Flux_pres_work_Cell_z=new double [N_Sam_CELLS];
	Q_heat_flux_x=new double[N_Sam_CELLS];
	Q_heat_flux_y=new double[N_Sam_CELLS];
	Q_heat_flux_z=new double[N_Sam_CELLS];
	Temp_Of_Sam_Cell=new double[N_Sam_CELLS];
	M_of_Atoms_Passed_Sam_Cell=new double[N_Sam_CELLS];
	N_of_Atoms_Passed_Sam_Cell=new int[N_Sam_CELLS];
	for (i=0;i<N_Sam_CELLS;i++) {
		KKK1=i%N_Sam_CELL[0];
		KKK2=(i%(N_Sam_CELL[0]*N_Sam_CELL[1]))/N_Sam_CELL[0];
		KKK3=i/(N_Sam_CELL[0]*N_Sam_CELL[1]);
		X_Of_Sam_Cell[i]=X_Low_E+KKK1*Dx_s[0]+Dx_s[0]/2;
		Y_Of_Sam_Cell[i]=Y_Low_E+KKK2*Dx_s[1]+Dx_s[1]/2;
		Z_Of_Sam_Cell[i]=Z_Low_E+KKK3*Dx_s[2]+Dx_s[2]/2;
		KE_Of_Sam_Cell[i]=0;
	
        Flux_KE_Of_X_Of_Sam_Cell[i]=0;
	    Flux_KE_Of_Y_Of_Sam_Cell[i]=0;
		Flux_KE_Of_Z_Of_Sam_Cell[i]=0;
	    Flux_pot_of_Cell_x[i]=0;
	    Flux_pot_of_Cell_y[i]=0;
		Flux_pot_of_Cell_z[i]=0;
	    Flux_pres_work_Cell_x[i]=0;
	    Flux_pres_work_Cell_y[i]=0;
		Flux_pres_work_Cell_z[i]=0;
	
		Vel_Of_Sam_Cell[i][0]=0;
		Vel_Of_Sam_Cell[i][1]=0;
		Vel_Of_Sam_Cell[i][2]=0;
		M_of_Atoms_Passed_Sam_Cell[i]=0;
		N_of_Atoms_Passed_Sam_Cell[i]=0;
		Density_Of_Sam_Cell[i]=0;
	}
	switch (Smp_INIT_MODE) {
		case 'I': {			// Initializing
		//	srand(B(0));
			} 
            break;
        case 'R': 	Read_Sam_Restart(); // read from arestart file.
             break;
        default : cout << "Error Input Mode!" << '\n';
    }
}
/***********************************************************************/
/***********************************************************************/
void Make_Sam_Cell_Member()
{
	int i1;
	for (i1=0;i1<N_Sam_CELLS;i1++) {
		Sam_CELL_MEMBER[i1][0]=0;
	}
	for (i1=0;i1<Num_Of_Atoms;i1++) {
		Sam_Cell_Det(i1) ;
		Sam_CELL_MEMBER[Atoms[i1].Sam_Cell][0]++;
		Sam_CELL_MEMBER[Atoms[i1].Sam_Cell][Sam_CELL_MEMBER[Atoms[i1].Sam_Cell][0]]=i1;
	}
}
/***********************************************************************/
/***********************************************************************/
void  Sam_Cell_Det(int i1) 
{
int 	KK=0,KK1,KK2,KK3;
		KK1=int(floor((Atoms[i1].R_Vec[0]-X_Low_E)/Dx_s[0]));
		if (KK1==-1) KK1=0;
		else if (KK1==N_Sam_CELL[0]) KK1=N_Sam_CELL[0]-1;
		KK2=int(floor((Atoms[i1].R_Vec[1]-Y_Low_E)/Dx_s[1]));
		if (KK2==-1) KK2=0;
		else if (KK2==N_Sam_CELL[1]) KK2=N_Sam_CELL[1]-1;
		
		KK3=int(floor((Atoms[i1].R_Vec[2]-Z_Low_E)/Dx_s[2]));
		if (KK3==-1) KK3=0;
		else if (KK3==N_Sam_CELL[2]) KK3=N_Sam_CELL[2]-1;
		Atoms[i1].Sam_Cell=N_Sam_CELL[0]*N_Sam_CELL[1]*KK3+N_Sam_CELL[0]*KK2+KK1;;
}
/***********************************************************************/
/***********************************************************************/
void Sampling() 
{
	int	i,i1,i3,NN;
	double Vel[3],KE_of_step_of_Cell,Flux_pot_of_step_of_Cell_x,Flux_pres_of_step_of_Cell_x,vx2,vy2,vz2,mm,MM,a_mass;
    double Flux_pot_of_step_of_Cell_y,Flux_pot_of_step_of_Cell_z,Flux_pres_of_step_of_Cell_y,Flux_pres_of_step_of_Cell_z,Flux_KE_of_step_of_Cell_x,Flux_KE_of_step_of_Cell_y,Flux_KE_of_step_of_Cell_z;
	double Tot_M_Passed;
	for (i=0;i<N_Sam_CELLS;i++) {
		NN=0;
		Vel[0]=Vel[1]=Vel[2]=KE_of_step_of_Cell=Flux_pot_of_step_of_Cell_y=Flux_pres_of_step_of_Cell_y=Flux_KE_of_step_of_Cell_x=Flux_KE_of_step_of_Cell_y=MM=0;
		Flux_pot_of_step_of_Cell_x=Flux_pres_of_step_of_Cell_x=Flux_pot_of_step_of_Cell_z=Flux_pres_of_step_of_Cell_z=Flux_KE_of_step_of_Cell_z=0;
		for (i1=1;i1<Sam_CELL_MEMBER[i][0]+1;i1++) {
			i3=Sam_CELL_MEMBER[i][i1];
			mm=AtomMass[Atoms[i3].tp-1];
			Vel[0]+=Atoms[i3].V_Vec[0]; 
			Vel[1]+=Atoms[i3].V_Vec[1];
			Vel[2]+=Atoms[i3].V_Vec[2];
			MM+=mm;
			NN++;
		}
		if (Sam_CELL_MEMBER[i][0]!=0) {
			Vel_Of_Sam_Cell[i][0]=(Vel_Of_Sam_Cell[i][0]*N_of_Atoms_Passed_Sam_Cell[i]+Vel[0]);
			Vel_Of_Sam_Cell[i][1]=(Vel_Of_Sam_Cell[i][1]*N_of_Atoms_Passed_Sam_Cell[i]+Vel[1]);
			Vel_Of_Sam_Cell[i][2]=(Vel_Of_Sam_Cell[i][2]*N_of_Atoms_Passed_Sam_Cell[i]+Vel[2]);
		}
		M_of_Atoms_Passed_Sam_Cell[i]+=MM;
		N_of_Atoms_Passed_Sam_Cell[i]+=NN;
		if (Sam_CELL_MEMBER[i][0]!=0) {
			Vel_Of_Sam_Cell[i][0]/=N_of_Atoms_Passed_Sam_Cell[i];
			Vel_Of_Sam_Cell[i][1]/=N_of_Atoms_Passed_Sam_Cell[i];
			Vel_Of_Sam_Cell[i][2]/=N_of_Atoms_Passed_Sam_Cell[i];
		}
		for (i1=1;i1<Sam_CELL_MEMBER[i][0]+1;i1++) {
			i3=Sam_CELL_MEMBER[i][i1];
			vx2=Atoms[i3].V_Vec[0]*Atoms[i3].V_Vec[0];
			vy2=Atoms[i3].V_Vec[1]*Atoms[i3].V_Vec[1];
			vz2=Atoms[i3].V_Vec[2]*Atoms[i3].V_Vec[2];

			a_mass=AtomMass[Atoms[i3].tp-1];
			KE_of_step_of_Cell+=0.5*a_mass*(vx2+vy2+vz2);
			
			Flux_KE_of_step_of_Cell_x+=(a_mass*(vx2+vy2+vz2))*Atoms[i3].V_Vec[0];
            Flux_KE_of_step_of_Cell_y+=(a_mass*(vx2+vy2+vz2))*Atoms[i3].V_Vec[1];
			Flux_KE_of_step_of_Cell_z+=(a_mass*(vx2+vy2+vz2))*Atoms[i3].V_Vec[2];
			Flux_pot_of_step_of_Cell_x+=Atoms[i3].Pot*Atoms[i3].V_Vec[0];
			Flux_pot_of_step_of_Cell_y+=Atoms[i3].Pot*Atoms[i3].V_Vec[1];
			Flux_pot_of_step_of_Cell_z+=Atoms[i3].Pot*Atoms[i3].V_Vec[2];
            Flux_pres_of_step_of_Cell_x+=Atoms[i3].Pres_Work*Atoms[i3].V_Vec[0];
			Flux_pres_of_step_of_Cell_y+=Atoms[i3].Pres_Work*Atoms[i3].V_Vec[1];
			Flux_pres_of_step_of_Cell_z+=Atoms[i3].Pres_Work*Atoms[i3].V_Vec[2];

		}
		if (Sam_CELL_MEMBER[i][0]!=0) {
			KE_Of_Sam_Cell[i]=(KE_Of_Sam_Cell[i]+KE_of_step_of_Cell);
			
			Flux_KE_Of_X_Of_Sam_Cell[i]+=Flux_KE_of_step_of_Cell_x;
            Flux_KE_Of_Y_Of_Sam_Cell[i]+=Flux_KE_of_step_of_Cell_y;
			Flux_KE_Of_Z_Of_Sam_Cell[i]+=Flux_KE_of_step_of_Cell_z;
            Flux_pot_of_Cell_x[i]+=Flux_pot_of_step_of_Cell_x;
			Flux_pot_of_Cell_y[i]+=Flux_pot_of_step_of_Cell_y;
			Flux_pot_of_Cell_z[i]+=Flux_pot_of_step_of_Cell_z;
            Flux_pres_work_Cell_x[i]+=Flux_pres_of_step_of_Cell_x;
			Flux_pres_work_Cell_y[i]+=Flux_pres_of_step_of_Cell_y;
			Flux_pres_work_Cell_z[i]+=Flux_pres_of_step_of_Cell_z;
		} 
	}
	Sam_Step++;
	if (Sam_Step%50==0) {
		Tot_M_Passed=0;
		for (i=0;i<N_Sam_CELLS;i++) {
			Tot_M_Passed+=M_of_Atoms_Passed_Sam_Cell[i];
		}
		for (i=0;i<N_Sam_CELLS;i++) {
			if(N_of_Atoms_Passed_Sam_Cell[i]!=0) {
				Temp_Of_Sam_Cell[i]=(KE_Of_Sam_Cell[i]-0.5*(M_of_Atoms_Passed_Sam_Cell[i])*(Vel_Of_Sam_Cell[i][0]*Vel_Of_Sam_Cell[i][0]+Vel_Of_Sam_Cell[i][1]*Vel_Of_Sam_Cell[i][1]+Vel_Of_Sam_Cell[i][2]*Vel_Of_Sam_Cell[i][2]))/(1.5*Boltzmann_const*N_of_Atoms_Passed_Sam_Cell[i]);
			
			    Q_heat_flux_x[i]=(1/(2*(Dx_s[0]*Dx_s[1]*Dx_s[2])*N_of_Atoms_Passed_Sam_Cell[i]))*(Flux_KE_Of_X_Of_Sam_Cell[i]+Flux_pot_of_Cell_x[i]-Flux_pres_work_Cell_x[i]);
			    Q_heat_flux_y[i]=(1/(2*(Dx_s[0]*Dx_s[1]*Dx_s[2])*N_of_Atoms_Passed_Sam_Cell[i]))*(Flux_KE_Of_Y_Of_Sam_Cell[i]+Flux_pot_of_Cell_y[i]-Flux_pres_work_Cell_y[i]);
				Q_heat_flux_z[i]=(1/(2*(Dx_s[0]*Dx_s[1]*Dx_s[2])*N_of_Atoms_Passed_Sam_Cell[i]))*(Flux_KE_Of_Z_Of_Sam_Cell[i]+Flux_pot_of_Cell_z[i]-Flux_pres_work_Cell_z[i]);
			
			} else {
				Temp_Of_Sam_Cell[i]=0;
				Q_heat_flux_x[i]=0;
				Q_heat_flux_y[i]=0;
				Q_heat_flux_z[i]=0;
			}
			Density_Of_Sam_Cell[i]=(M_of_Atoms_Passed_Sam_Cell[i])/Tot_M_Passed*M_system/Volume_Of_Sam_Cell;

		}
		Write_Sampled_Data();
	} 
}
/***********************************************************************/
/***********************************************************************/
void Write_Sampled_Data()
{
	int i;
	ofstream Sam_Wrt("Sampled_Data.plt");
	Sam_Wrt << "VARIABLES=" <<'"' << "X" << '"' <<"," << '"'<< "Y"<< '"' <<"," << '"'<< "Z" << '"' << "," << '"' << "U" << '"'<< "," <<'"' << "V"<< '"'<< "," <<'"' << "W"<<'"' <<"," <<'"' << "Rho"<<'"' <<"," <<'"' << "T"<<'"' <<"," <<'"' << "Qx"<<'"' <<"," <<'"' <<"Qy"<<'"' <<"," <<'"' <<"Qz"<<'"' <<'\n';
	Sam_Wrt << "ZONE  I=" <<N_Sam_CELL[0] <<" J="<<N_Sam_CELL[1]<<" K="<<N_Sam_CELL[2] << " F=" << '"' << "POINT" << '"' << '\n';
	for (i=0;i<N_Sam_CELLS;i++) {
		Sam_Wrt << X_Of_Sam_Cell[i] << "		" << Y_Of_Sam_Cell[i]<< "		" << Z_Of_Sam_Cell[i] << "		" << Vel_Of_Sam_Cell[i][0] << "		" << Vel_Of_Sam_Cell[i][1]<< "		" << Vel_Of_Sam_Cell[i][2] <<  "	" << Density_Of_Sam_Cell[i] <<  "    	"  << Temp_Of_Sam_Cell[i]<<  "    	"<< Q_heat_flux_x[i]<< "       "<<Q_heat_flux_y[i]<< "       "<<Q_heat_flux_z[i]<<'\n';
	}
	Sam_Wrt.close();
	if (Sam_Step%5000==0) {
		Write_Sam_Restart();
	}
}
/***********************************************************************/
/***********************************************************************/
void Write_Sam_Restart()
{
	int i;
	ofstream Sam_rst("Sampled_Restart.dat");
	Sam_rst << "  M_tot			N_tot		K		U		V	 W		F_K_X	F_K_Y	F_K_Z	F_P_X	F_P_Y	F_P_Z	F_W_X	F_W_Y	F_W_Z" <<'\n';
	for (i=0;i<N_Sam_CELLS;i++) {
		Sam_rst	<< i <<"    " << M_of_Atoms_Passed_Sam_Cell[i] << "   "<<	N_of_Atoms_Passed_Sam_Cell[i]  << "   "<< KE_Of_Sam_Cell[i]  << "   "<< Vel_Of_Sam_Cell[i][0]	<< "   "<<	Vel_Of_Sam_Cell[i][1]<< "   "<<	Vel_Of_Sam_Cell[i][2];
		Sam_rst	<< "   "<<	Flux_KE_Of_X_Of_Sam_Cell[i]<< "   "<<	Flux_KE_Of_Y_Of_Sam_Cell[i]<<  "    "<<	Flux_KE_Of_Z_Of_Sam_Cell[i]<<  "    "<<	Flux_pot_of_Cell_x[i]<<  "    "  <<	Flux_pot_of_Cell_y[i]<<  "    "  <<	Flux_pot_of_Cell_z[i]<< "   "<<	Flux_pres_work_Cell_x[i] <<  "    " <<	Flux_pres_work_Cell_y[i] <<  "    " <<	Flux_pres_work_Cell_z[i]<< '\n';
	}
	Sam_rst << "complete" ;
	Sam_rst.close();
	ofstream Sam_rst2("Sampled_Restart_Copy.dat");
	Sam_rst2 << "  M_tot			N_tot		K		U		V	 W		F_K_X	F_K_Y	F_K_Z	F_P_X	F_P_Y	F_P_Z	F_W_X	F_W_Y	F_W_Z" <<'\n';
	for (i=0;i<N_Sam_CELLS;i++) {
		Sam_rst2	<< i <<"    " << M_of_Atoms_Passed_Sam_Cell[i] << "   "<<	N_of_Atoms_Passed_Sam_Cell[i]  << "   "<< KE_Of_Sam_Cell[i]  << "   "<< Vel_Of_Sam_Cell[i][0]	<< "   "<<	Vel_Of_Sam_Cell[i][1]<< "   "<<	Vel_Of_Sam_Cell[i][2];
		Sam_rst2	<< "   "<<	Flux_KE_Of_X_Of_Sam_Cell[i]<< "   "<<	Flux_KE_Of_Y_Of_Sam_Cell[i]<<  "    "<<	Flux_KE_Of_Z_Of_Sam_Cell[i]<<  "    "<<	Flux_pot_of_Cell_x[i]<<  "    "  <<	Flux_pot_of_Cell_y[i]<<  "    "  <<	Flux_pot_of_Cell_z[i]<< "   "<<	Flux_pres_work_Cell_x[i] <<  "    " <<	Flux_pres_work_Cell_y[i] <<  "    " <<	Flux_pres_work_Cell_z[i]<< '\n';
	}
	Sam_rst2 << "complete" ;
	Sam_rst2.close();
}


/***********************************************************************/
/***********************************************************************/
void	Immobalize_System()
{
	int	i;
	double M_sys=0,a_mass;
	v_mean[0]=0;
	v_mean[1]=0;
	v_mean[2]=0;
	for (i=0;i<Num_Of_Atoms;i++) {
		a_mass=AtomMass[Atoms[i].tp-1];
		v_mean[0]+=(a_mass*Atoms[i].V_Vec[0]);
		v_mean[1]+=(a_mass*Atoms[i].V_Vec[1]);
		v_mean[2]+=(a_mass*Atoms[i].V_Vec[2]);
		M_sys+=a_mass;
	}
	v_mean[0]=v_mean[0]/M_sys;
	v_mean[1]=v_mean[1]/M_sys;
	v_mean[2]=v_mean[2]/M_sys;
	for (i=0;i<Num_Of_Atoms;i++) {
		Atoms[i].V_Vec[0]-=v_mean[0];
		Atoms[i].V_Vec[1]-=v_mean[1];
		Atoms[i].V_Vec[2]-=v_mean[2];
	}
}
/***********************************************************************/
/***********************************************************************/
void	Write_Trajectory_Restart()
{
	ofstream rst_1("Restart_1.dat");
	int i;
	rst_1 << "Time"  << "		" << Time << '\n'; 
	rst_1 << "Xi"  << "		" << Xi << '\n'; 
	rst_1 << "Xip"<< "		"  << Xip << '\n';
	rst_1 << "Step" << "      "  << StepNum << '\n';
	rst_1 << "Num_of_atoms" << "      "  << Num_Of_Atoms << '\n';
	if ((Integ_Method=='V') ||(Integ_Method=='v') ) {
		rst_1 << "  n   " << "   X    "  <<  "   Y   " <<  "   Z  " <<  "  u  "    <<  "  v  " <<  "  w  "    << " ax   "  << "  ay  " << "  az  " << "  tp  "<< '\n';
    	for (i=0;i<Num_Of_Atoms;i++) {
			rst_1 << i <<"  " <<  Atoms[i].R_Vec[0]<< "	" << Atoms[i].R_Vec[1]<< "	" << Atoms[i].R_Vec[2] <<  "	" <<Atoms[i].V_Vec[0] <<  "	" << Atoms[i].V_Vec[1] <<  "	"<< Atoms[i].V_Vec[2] <<  "	" << Atoms[i].A_Vec[0] <<  "	" << Atoms[i].A_Vec[1]<<  "	" << Atoms[i].A_Vec[2]<<  "	" << Atoms[i].tp<<  '\n';
		}	
	} 
	rst_1 << "Complete";
	rst_1.close();
	ofstream rst_2("Restart_2.dat");
	rst_2 << "Time"  << "		" << Time << '\n'; 
	rst_2 << "Xi"  << "		" << Xi << '\n'; 
	rst_2 << "Xip"<< "		"  << Xip << '\n';
	rst_2 << "Step" << "      "  << StepNum << '\n';
	rst_1 << "Num_of_atoms" << "      "  << Num_Of_Atoms << '\n';
	if ((Integ_Method=='V') ||(Integ_Method=='v') ) {
		rst_2 << "  n   " << "   X    "  <<  "   Y   " <<  "   Z  " <<  "  u  "    <<  "  v  " <<  "  w  "    << " ax   "  << "  ay  " << "  az  " <<  "  tp  " <<'\n';
    	for (i=0;i<Num_Of_Atoms;i++) {
			rst_2 << i <<"  " <<  Atoms[i].R_Vec[0]<< "	" << Atoms[i].R_Vec[1]<< "	" << Atoms[i].R_Vec[2] <<  "	" <<Atoms[i].V_Vec[0] <<  "	" << Atoms[i].V_Vec[1] <<  "	"<< Atoms[i].V_Vec[2] <<  "	" << Atoms[i].A_Vec[0] <<  "	" << Atoms[i].A_Vec[1]<<  "	" << Atoms[i].A_Vec[2]<<  "	" << Atoms[i].tp<<  '\n';
		}	
	} 
	rst_2 << "Complete";
	rst_2.close();
}
/***********************************************************************/
/***********************************************************************/
void Read__Trajectory_Restart() 
{
	int i,j;
	char D[40];
	ifstream Inp_1("Restart_1.dat");
	Inp_1 >> D  >> Time;
	int hh=2;
	Inp_1 >> D  >> Xi;
	Inp_1 >> D  >> Xip;
	Inp_1 >> D  >> StepNum;
	Inp_1 >> D  >> Num_Of_Atoms;
	if ((Integ_Method=='V') ||(Integ_Method=='v') ) {
		for (i=0;i<11;i++) Inp_1 >> D;
		for(i=0;i<Num_Of_Atoms;i++) {
			Inp_1 >> j >>   Atoms[i].R_Vec[0] >>  Atoms[i].R_Vec[1]  >> Atoms[i].R_Vec[2] >> Atoms[i].V_Vec[0]  >>  Atoms[i].V_Vec[1]  >> Atoms[i].V_Vec[2] >>  Atoms[i].A_Vec[0]  >>  Atoms[i].A_Vec[1]  >> Atoms[i].A_Vec[2] >>  Atoms[i].tp ;
		}
	} 
	Inp_1.close();
}
/***********************************************************************/
/***********************************************************************/
void Read_Sam_Restart() 
{
	int i,j;
	char D[40];
	ifstream Sam_rst("Sampled_Restart.dat");
	Sam_rst >>  D  	>> D  >> D  >> D >> D  >> D >> D  >> D   >> D   >> D   >> D >>D >> D >> D >>D ;
	for (i=0;i<N_Sam_CELLS;i++) {
		Sam_rst	>> j >> M_of_Atoms_Passed_Sam_Cell[i] >> N_of_Atoms_Passed_Sam_Cell[i] >> KE_Of_Sam_Cell[i] >> Vel_Of_Sam_Cell[i][0]	>> Vel_Of_Sam_Cell[i][1] >> Vel_Of_Sam_Cell[i][2];
		Sam_rst	>> Flux_KE_Of_X_Of_Sam_Cell[i] >> 	Flux_KE_Of_Y_Of_Sam_Cell[i] >>	Flux_KE_Of_Z_Of_Sam_Cell[i] >>	Flux_pot_of_Cell_x[i]  >>	Flux_pot_of_Cell_y[i] >>	Flux_pot_of_Cell_z[i]>>  Flux_pres_work_Cell_x[i]   >> Flux_pres_work_Cell_y[i]>> Flux_pres_work_Cell_z[i] ;
	}
	Sam_rst.close();
}
/***********************************************************************/
/***********************************************************************/
void Propeller_Force()
{
	double  *Scale_Factor, Tot_Drag = 0.0;
	int I, X, Y;

	Scale_Factor = new double[Num_Of_Thrust_Sets];

	for(int i=0; i<Num_Of_Thrust_Sets; i++)
	{
		X = Thrust_Set[i].Atom_No_Direction_1 - 1;
		Y = Thrust_Set[i].Atom_No_Direction_2 - 1;

		Thrust_Set[i].xx = Atoms[X].R_Vec[0] - Atoms[X].R_Vec[0];
		Thrust_Set[i].yy = Atoms[X].R_Vec[1] - Atoms[X].R_Vec[1];
		Thrust_Set[i].zz = Atoms[X].R_Vec[2] - Atoms[X].R_Vec[2];

		Scale_Factor[i] = (Thrust_Set[i].F) / (sqrt(Thrust_Set[i].xx * Thrust_Set[i].xx + Thrust_Set[i].yy * Thrust_Set[i].yy + Thrust_Set[i].zz * Thrust_Set[i].zz));

		Thrust_Set[i].fx = Scale_Factor[i] * ((Thrust_Set[i].xx));
		Thrust_Set[i].fy = Scale_Factor[i] * ((Thrust_Set[i].yy));
		Thrust_Set[i].fz = Scale_Factor[i] * ((Thrust_Set[i].zz));

		I = Thrust_Set[i].Atom_No - 1;

		Atoms[I].force[0] += Thrust_Set[i].fx;
		Atoms[I].force[1] += Thrust_Set[i].fy;
		Atoms[I].force[2] += Thrust_Set[i].fz;
	}
}
/***********************************************************************/
/***********************************************************************/
void Ref_Del_List () 
{
	int i;
 	for (i=0;i<Num_Of_Atoms;i++) Atoms[i].Delete_List=0;
}
/***********************************************************************/
/***********************************************************************/

void Add_Atom(double xx,double yy,double zz,double uu,double vv,double ww,int tp)
{
	double MM;
	MM=AtomMass[Atoms[Num_Of_Atoms-1].tp-1];

	Atoms[Num_Of_Atoms-1].R_Vec[0]=xx;
	Atoms[Num_Of_Atoms-1].R_Vec[1]=yy;
	Atoms[Num_Of_Atoms-1].R_Vec[2]=zz;

	Atoms[Num_Of_Atoms-1].V_Vec[0]=uu;
	Atoms[Num_Of_Atoms-1].V_Vec[1]=vv;
	Atoms[Num_Of_Atoms-1].V_Vec[2]=ww;
	Atoms[Num_Of_Atoms-1].tp=tp;
	M_system+=MM;

	Create_Neighbourhood_of_Atom(Num_Of_Atoms-1);
	Atoms[Num_Of_Atoms-1].force[0]=0;
	Atoms[Num_Of_Atoms-1].force[1]=0;
	Atoms[Num_Of_Atoms-1].force[2]=0;

	LJ_COUL_of_Atom(Num_Of_Atoms-1);

	Atoms[Num_Of_Atoms-1].A_Vec[0]=Atoms[Num_Of_Atoms-1].force[0]/MM;
	Atoms[Num_Of_Atoms-1].A_Vec[1]=Atoms[Num_Of_Atoms-1].force[1]/MM;
	Atoms[Num_Of_Atoms-1].A_Vec[2]=Atoms[Num_Of_Atoms-1].force[2]/MM;
}
/***********************************************************************/
/***********************************************************************/
void Insert_Atom_in_Cell(int I, int N,double xx1,double xx2,double yy1,double yy2,double zz1,double zz2,double velx,double vely,double velz, double temp,int tp)
{
	int i;
	for (i=0;i<N;i++) {
		Usher_in_Cell(I,xx1,xx2,yy1,yy2,zz1,zz2,velx,vely,velz,temp,tp);
	}
}
/***********************************************************************/
/***********************************************************************/
void Usher_in_Cell(int I,double xx1,double xx2,double yy1,double yy2,double zz1,double zz2,double velx,double vely,double velz, double temp,int tp)
{
	double x_ij,y_ij,z_ij,r_ij2,r_ij2_inv,r_ij6_inv,r_ij12_inv,Force,r_ij,x_ij2,y_ij2,z_ij2;
	int KK1,KK2,KK3,i,jj,j;
	double X_L,X_H,Y_L,Y_H,Z_L,Z_H,xx,yy,zz,uu,vv,ww,U0,du,F,U_ovlp,Ds,Ds1,Ds_ovlp;
	double sigma,epsilon;
	int NNN,Reinit;
	double dx,dy,dz,RHO;
	double s_T;
	int numofNeighbors;
	int N_try=20;
	int Ins_OK=1;
	RHO=CELL_MEMBER[I][0]*sigma_max*sigma_max*sigma_max/(Dx[0]*Dx[1]*Dx[2]);
	if (RHO==0) RHO=0.001;
	Ds1=0.1*pow(RHO,(-3.0/2));
	if (Ds1>Dx[0]/4)  Ds1=Dx[0]/4;
	KK1=I%N_CELL[0];
	KK2=(I%(N_CELL[0]*N_CELL[1]))/N_CELL[0];
	KK3=I/(N_CELL[0]*N_CELL[1]);
	X_L=X_Low_E+KK1*Dx[0];
	Y_L=Y_Low_E+KK2*Dx[1];
	Z_L=Z_Low_E+KK3*Dx[2];
	X_H=X_L+Dx[0];
	Y_H=Y_L+Dx[1];
	Z_H=Z_L+Dx[2];
	xx=Uniform_Random(X_L,X_H);
	yy=Uniform_Random(Y_L,Y_H);
	zz=Uniform_Random(Z_L,Z_H);
	Num_Of_Atoms++;
	Atoms[Num_Of_Atoms-1].tp=tp;
	if (CELL_MEMBER[I][0]<0) {
		cout << "ERRRRRRRRRRRRRRRRRRor"<<'\n';
	}
	CELL_MEMBER[I][0]++;
	CELL_MEMBER[I][CELL_MEMBER[I][0]]=Num_Of_Atoms-1;
	du=10.0;
	U0=-0.0005;
	NNN=0;
	double 	Res=0.0001;
	int 	NN=0;
	U_ovlp=1.0;
	int Impossility_Fac=10000;
	Ds_ovlp=0.9*Sigma_MAT[tp-1][tp-1]-pow((4*Epsilon_MAT[tp-1][tp-1]/U_ovlp),(1.0/12))*Sigma_MAT[tp-1][tp-1];
	while (du>Res) {
		U0=U0+0.0001;
		NN=0;
		while (NN<N_try) {
			Atoms[Num_Of_Atoms-1].R_Vec[0]=xx;
			Atoms[Num_Of_Atoms-1].R_Vec[1]=yy;
			Atoms[Num_Of_Atoms-1].R_Vec[2]=zz;
			Create_Neighbourhood_of_Atom(Num_Of_Atoms-1);
			Atoms[Num_Of_Atoms-1].force[0]=0;
			Atoms[Num_Of_Atoms-1].force[1]=0;
			Atoms[Num_Of_Atoms-1].force[2]=0;
			i=Num_Of_Atoms-1;
			numofNeighbors=Atoms[i].NeighborList[0];
			Atoms[i].Pot=0;
			for (jj=1;jj<numofNeighbors+1;jj++) {
				j=Atoms[i].NeighborList[jj];
				z_ij=(Atoms[i].R_Vec[2]-Atoms[j].R_Vec[2]);
				y_ij=(Atoms[i].R_Vec[1]-Atoms[j].R_Vec[1]);
				x_ij=(Atoms[i].R_Vec[0]-Atoms[j].R_Vec[0]);
				if ((Boundary_Type.x=='p') |(Boundary_Type.x=='P')) {
					if ((x_ij)>X_H_L/2)  x_ij=x_ij-X_H_L;
					else if ((x_ij)<-X_H_L/2) x_ij=x_ij+X_H_L;
				}
				if ((Boundary_Type.y=='p') |(Boundary_Type.y=='P')) {
					if ((y_ij)>Y_H_L/2) y_ij=y_ij-Y_H_L;
					else if ((y_ij)<-Y_H_L/2) y_ij=y_ij+Y_H_L;
				}
				if ((Boundary_Type.z=='p') |(Boundary_Type.z=='P')) {
					if ((z_ij)>Z_H_L/2) z_ij=z_ij-Z_H_L;
					else if ((z_ij)<-Z_H_L/2) z_ij=z_ij+Z_H_L;
				}
				x_ij2=x_ij*x_ij;
				y_ij2=y_ij*y_ij;
				z_ij2=z_ij*z_ij;
				r_ij2=(x_ij2+y_ij2+z_ij2);
				r_ij=sqrt(r_ij2);
				sigma=Sigma_MAT[Atoms[i].tp-1][Atoms[j].tp-1];
				R_Cut_Off=R_Cut_Off_R*sigma;
				if ((r_ij<R_Cut_Off)) { 
					r_ij2_inv=sigma*sigma/r_ij2;
					r_ij6_inv=r_ij2_inv*r_ij2_inv*r_ij2_inv;
					r_ij12_inv=r_ij6_inv*r_ij6_inv;
					r_ij2_inv=r_ij2_inv/(sigma*sigma);
					epsilon=Epsilon_MAT[Atoms[i].tp-1][Atoms[j].tp-1];
					Force=48*epsilon*r_ij2_inv*(r_ij12_inv-0.5*r_ij6_inv);
					Atoms[i].force[0]+=x_ij*Force;
					Atoms[i].force[1]+=y_ij*Force;
					Atoms[i].force[2]+=z_ij*Force;
					Atoms[i].Pot+=4*epsilon*(r_ij12_inv-r_ij6_inv);
			   }
			}
			du=fabs((Atoms[i].Pot-U0)/U0);
			F=sqrt(Atoms[i].force[0]*Atoms[i].force[0]+Atoms[i].force[1]*Atoms[i].force[1]+Atoms[i].force[2]*Atoms[i].force[2]);			
			Ds=Ds1;
			if (Atoms[i].Pot>U_ovlp) {
				Ds_ovlp=0.9*Sigma_MAT[tp-1][tp-1]-pow((4*Epsilon_MAT[tp-1][tp-1]/Atoms[i].Pot),(1.0/12))*Sigma_MAT[tp-1][tp-1];
				Ds=Ds_ovlp;
			}else {
				if (((Atoms[i].Pot-U0)/F)<Ds1) {
					Ds=(Atoms[i].Pot-U0)/F;
				}
			}
			if (F<.000001) du=0.0;
			else {
				dx=Atoms[i].force[0]/F*Ds;
				dy=Atoms[i].force[1]/F*Ds;
				dz=Atoms[i].force[2]/F*Ds;
				xx+=dx;
				yy+=dy;
				zz+=dz;
			}			
			NNN++;
			Reinit=0;
			if (xx>X_H)  Reinit=1;
			else if (xx<X_L)  Reinit=1;
			else if (yy>Y_H)  Reinit=1;
			else if (yy<Y_L)  Reinit=1;
			else if (zz>Z_H)  Reinit=1;
			else if (zz<Z_L)  Reinit=1;

			if (NNN>200)	  Reinit=1;						
			if  (Reinit==1) {
				NNN=0;
				xx=Uniform_Random(X_L,X_H);
				yy=Uniform_Random(Y_L,Y_H);
				zz=Uniform_Random(Z_L,Z_H);
			}
			NN++;
			if (du<Res) break;
		}
		if(U0>U_ovlp/Impossility_Fac)  {
			Ins_OK=0;
			break;
		}			
	}
	if (Ins_OK==1) {
		s_T=sqrt(Boltzmann_const*temp/AtomMass[tp-1]);
		uu=Normal_Random(0,s_T);
		vv=Normal_Random(0,s_T);
		ww=Normal_Random(0,s_T);	

		if ((xx>xx1) && (xx<xx2)) {
			if ((yy>yy1) && (yy<yy2)) {
				if ((zz>zz1) && (zz<zz2)) {
					uu+=velx;
					vv+=vely;
					ww+=velz;
				}
			}
		}
		Add_Atom(xx,yy,zz,uu,vv,ww,tp);	
	} else {
		Num_Of_Atoms--;
		CELL_MEMBER[I][0]--;		
	}
}

/***********************************************************************/
/***********************************************************************/
void Delete_Atoms() 
{
	int nmn=(2*MMM+1)*(2*MMM+1)*(2*MMM+1);
	int i,i2,k3,k4,i6,i1,i5;
	int Mx_N=Num_Of_Atoms;
	
	/*************** DELETING BONDS & ANGLES ***************/
	
	int atm1,atm2,atm3;
	for (i=0;i<Num_Of_Bonds;i++) {
		atm1=bondlist[i][0];
		atm2=bondlist[i][1];
		if ((Atoms[atm1].Delete_List==1)|(Atoms[atm2].Delete_List==1)){
			Atoms[atm1].Delete_List=1;
			Atoms[atm2].Delete_List=1;
			bondlist[i][0]=bondlist[Num_Of_Bonds-1][0];
			bondlist[i][1]=bondlist[Num_Of_Bonds-1][1];
			bondlist[i][2]=bondlist[Num_Of_Bonds-1][2];
			Num_Of_Bonds--;
		}
	}

	for (i=0;i<Num_Of_Angles;i++) {
		atm1=anglelist[i][0];
		atm2=anglelist[i][1];
		atm3=anglelist[i][2];
		if (((Atoms[atm1].Delete_List==1)|(Atoms[atm2].Delete_List==1))|(Atoms[atm3].Delete_List==1)){
			Atoms[atm1].Delete_List=1;
			Atoms[atm2].Delete_List=1;
			Atoms[atm3].Delete_List=1;
			anglelist[i][0]=anglelist[Num_Of_Angles-1][0];
			anglelist[i][1]=anglelist[Num_Of_Angles-1][1];
			anglelist[i][2]=anglelist[Num_Of_Angles-1][2];
			anglelist[i][2]=anglelist[Num_Of_Angles-1][2];
			Num_Of_Angles--;
		}
	}
	/****************************************************/
	i=0;
	while (i<Num_Of_Atoms-1) {
		if (Atoms[i].Delete_List==1) {
			while (Atoms[i].Delete_List==1) {
				if (i!=Num_Of_Atoms-1) {
					Atoms[i].R_Vec[0]=Atoms[Num_Of_Atoms-1].R_Vec[0];
					Atoms[i].R_Vec[1]=Atoms[Num_Of_Atoms-1].R_Vec[1];
					Atoms[i].R_Vec[2]=Atoms[Num_Of_Atoms-1].R_Vec[2];
					Atoms[i].V_Vec[0]=Atoms[Num_Of_Atoms-1].V_Vec[0];
					Atoms[i].V_Vec[1]=Atoms[Num_Of_Atoms-1].V_Vec[1];
					Atoms[i].V_Vec[2]=Atoms[Num_Of_Atoms-1].V_Vec[2];
					Atoms[i].A_Vec[0]=Atoms[Num_Of_Atoms-1].A_Vec[0];
					Atoms[i].A_Vec[1]=Atoms[Num_Of_Atoms-1].A_Vec[1];
					Atoms[i].A_Vec[2]=Atoms[Num_Of_Atoms-1].A_Vec[2];
					Atoms[i].tp=Atoms[Num_Of_Atoms-1].tp;
					Atoms[i].Pres_Work=Atoms[Num_Of_Atoms-1].Pres_Work;
					Atoms[i].Inc_vec=Atoms[Num_Of_Atoms-1].Inc_vec;
					M_system-=AtomMass[Atoms[i].tp-1];
					i5=Atoms[i].Cell;
					for (k3=0;k3<nmn;k3++) {
						k4=Neighb_Cell[i5][k3];
						for (i6=1;i6<CELL_MEMBER[k4][0]+1;i6++) {
							i1=CELL_MEMBER[k4][i6];
							for (i2=1;i2<Atoms[i1].NeighborList[0]+1;i2++) {
								if (Atoms[i1].NeighborList[i2]==i) {
									Atoms[i1].NeighborList[i2]=Atoms[i1].NeighborList[Atoms[i1].NeighborList[0]];
									Atoms[i1].NeighborList[0]+=-1;
									break;
								}
							}
						}
					}
					for (i6=1;i6<CELL_MEMBER[i5][0]+1;i6++) {
						if (CELL_MEMBER[i5][i6]==i) {
							CELL_MEMBER[i5][i6]=CELL_MEMBER[i5][CELL_MEMBER[i5][0]];
							CELL_MEMBER[i5][0]+=-1;
							break;
						}
					}
					i5=Atoms[Num_Of_Atoms-1].Cell;
					for (k3=0;k3<nmn;k3++) {
						k4=Neighb_Cell[i5][k3];
						for (i6=1;i6<CELL_MEMBER[k4][0]+1;i6++) {
							i1=CELL_MEMBER[k4][i6];
							for (i2=1;i2<Atoms[i1].NeighborList[0]+1;i2++) {
								if (Atoms[i1].NeighborList[i2]==Num_Of_Atoms-1) {
									Atoms[i1].NeighborList[i2]=i;
									break;
								}
							}
						}
					}
					for (i6=1;i6<CELL_MEMBER[i5][0]+1;i6++) {
						if (CELL_MEMBER[i5][i6]==Num_Of_Atoms-1) {
							CELL_MEMBER[i5][i6]=i;
							break;
						}
					}
					int II=Atoms[Num_Of_Atoms-1].NeighborList[0]+1;
					for (i6=1;i6<II;i6++) {
						if (Atoms[Num_Of_Atoms-1].NeighborList[i6]==i) {
							II=i6;
							break;
						}
					}					
					for (i6=0;i6<II;i6++) {
						Atoms[i].NeighborList[i6]=Atoms[Num_Of_Atoms-1].NeighborList[i6];
					}
					if (II!=Atoms[Num_Of_Atoms-1].NeighborList[0]+1) {
						Atoms[i].NeighborList[0]=Atoms[Num_Of_Atoms-1].NeighborList[0]-1;
						for (i6=II;i6<Atoms[Num_Of_Atoms-1].NeighborList[0];i6++) {
							Atoms[i].NeighborList[i6]=Atoms[Num_Of_Atoms-1].NeighborList[i6+1];
						}
					}			
					Atoms[i].Delete_List=Atoms[Num_Of_Atoms-1].Delete_List;
					Atoms[i].Cell=Atoms[Num_Of_Atoms-1].Cell;
					Num_Of_Atoms--;
				} else {
					break;
				}
			}
		}
		i++;
	}
	i=Num_Of_Atoms-1;
	if (Atoms[i].Delete_List==1) {
		M_system-=AtomMass[Atoms[i].tp-1];
		i5=Atoms[i].Cell;
		for (k3=0;k3<nmn;k3++) {
			k4=Neighb_Cell[i5][k3];
			for (i6=1;i6<CELL_MEMBER[k4][0]+1;i6++) {
				i1=CELL_MEMBER[k4][i6];
				for (i2=1;i2<Atoms[i1].NeighborList[0]+1;i2++) {
					if (Atoms[i1].NeighborList[i2]==i) {
						Atoms[i1].NeighborList[i2]=Atoms[i1].NeighborList[Atoms[i1].NeighborList[0]];
						Atoms[i1].NeighborList[0]+=-1;
						break;
					}
				}
			}
		}
		for (i6=1;i6<CELL_MEMBER[i5][0]+1;i6++) {
			if (CELL_MEMBER[i5][i6]==i) {
				CELL_MEMBER[i5][i6]=CELL_MEMBER[i5][CELL_MEMBER[i5][0]];
				CELL_MEMBER[i5][0]+=-1;
				break;
			}
		}
		Atoms[i].Delete_List=0;
		Num_Of_Atoms--;
	}
}
/***********************************************************************/
/***********************************************************************/
void Inlets()
{
	int N1,i,j,k,m,ni;
	double Rho,temp;
	double Des_Rho;
	double velx,vely,velz,xx1,yy1,zz1,xx2,yy2,zz2;
	int tp;
	double s_T;
	double U_th,V_th,W_th,Kin_tot_Cell;
	double Fact,U_mean,V_mean,W_mean;
	double N_Ins;
	int *Cell_Inlet_Atoms=NULL;
	Cell_Inlet_Atoms=new int[Max_N_Of_Cell];
	for (ni=0;ni<Num_of_Inlets;ni++) {
		tp=Inlet[ni].atom_tp;
		Des_Rho=Inlet[ni].Rho;
		temp=Inlet[ni].temp;
		s_T=sqrt(Boltzmann_const*temp/AtomMass[tp-1]);
		velx=Inlet[ni].Vx;
		vely=Inlet[ni].Vy;
		velz=Inlet[ni].Vz;
		xx1=Inlet[ni].x1;
		xx2=Inlet[ni].x2;
		yy1=Inlet[ni].y1;
		yy2=Inlet[ni].y2;
		zz1=Inlet[ni].z1;
		zz2=Inlet[ni].z2;
		for (j=0;j<Inlet[ni].N_C;j++) {
			i=Inlet[ni].Cells[j];
			Cell_Inlet_Atoms[0]=0;
			U_mean=0;
			V_mean=0;
			W_mean=0;
			for (k=1;k<CELL_MEMBER[i][0]+1;k++) {
				m=CELL_MEMBER[i][k]; 
				if ((Atoms[m].R_Vec[0]>xx1) && (Atoms[m].R_Vec[0]<xx2)) {
					if ((Atoms[m].R_Vec[1]>yy1) && (Atoms[m].R_Vec[1]<yy2)) {
						if ((Atoms[m].R_Vec[2]>zz1) && (Atoms[m].R_Vec[2]<zz2)) {
							Cell_Inlet_Atoms[0]++;
							Cell_Inlet_Atoms[Cell_Inlet_Atoms[0]]=m;
							U_mean+=Atoms[m].V_Vec[0];
							V_mean+=Atoms[m].V_Vec[1];
							W_mean+=Atoms[m].V_Vec[2];
						}
					}
				}
			}
			if (Cell_Inlet_Atoms[0]>0) {
				U_mean=U_mean/Cell_Inlet_Atoms[0];
				V_mean=V_mean/Cell_Inlet_Atoms[0];
				W_mean=W_mean/Cell_Inlet_Atoms[0];
			}


			Kin_tot_Cell=0.0;
			for (k=1;k<Cell_Inlet_Atoms[0]+1;k++) {
				m=Cell_Inlet_Atoms[k];
				U_th=Atoms[m].V_Vec[0]-U_mean;
				V_th=Atoms[m].V_Vec[1]-V_mean;
				W_th=Atoms[m].V_Vec[2]-W_mean;
				Kin_tot_Cell+=AtomMass[Atoms[m].tp-1]*(U_th*U_th+V_th*V_th+W_th*W_th);
			}
			Kin_tot_Cell=Kin_tot_Cell/2.0;
			Fact=sqrt((1.5*Cell_Inlet_Atoms[0]*Boltzmann_const*temp)/Kin_tot_Cell);
			if (Cell_Inlet_Atoms[0]==1) {
				m=Cell_Inlet_Atoms[1];
				s_T=sqrt(Boltzmann_const*temp/AtomMass[Atoms[m].tp-1]);
				Atoms[m].V_Vec[0]=Normal_Random(velx,s_T);
				Atoms[m].V_Vec[1]=Normal_Random(vely,s_T);
				Atoms[m].V_Vec[2]=Normal_Random(velz,s_T);

			}else  {
				for (k=1;k<Cell_Inlet_Atoms[0]+1;k++) {
					m=Cell_Inlet_Atoms[k];
					U_th=Atoms[m].V_Vec[0]-U_mean;
					V_th=Atoms[m].V_Vec[1]-V_mean;
					W_th=Atoms[m].V_Vec[2]-W_mean;
					Atoms[m].V_Vec[0]=velx+U_th*Fact;
					Atoms[m].V_Vec[1]=vely+V_th*Fact;
					Atoms[m].V_Vec[2]=velz+W_th*Fact;
				}
				
			}
		}

		for (j=0;j<Inlet[ni].N_C;j++) {
			i=Inlet[ni].Cells[j];
			tp=Inlet[ni].atom_tp;
			Rho=CELL_MEMBER[i][0]*(sigma_max*sigma_max*sigma_max)/(Dx[0]*Dx[1]*Dx[2]);
			if (Rho<Des_Rho) {
				N_Ins=(Des_Rho-Rho)*(Dx[0]*Dx[1]*Dx[2])/(sigma_max*sigma_max*sigma_max)+Inlet[ni].Cell_N_Res[j];
				N1=int (floor(N_Ins));
				Inlet[ni].Cell_N_Res[j]=N_Ins-floor(N_Ins);
				Insert_Atom_in_Cell(i,N1,xx1,xx2,yy1,yy2,zz1,zz2,velx,vely,velz,temp,tp);
			}
		}
	}
	delete [] Cell_Inlet_Atoms;
	Cell_Inlet_Atoms=NULL;
}
/***********************************************************************/
/***********************************************************************/
int Belong_Cell_to_Inlet(int i,int i1)  
{
	double x_L,y_L,z_L,x_H,y_H,z_H;
	int KK1,KK2,KK3;
	KK1=i%N_CELL[0];
	KK2=(i%(N_CELL[0]*N_CELL[1]))/N_CELL[0];
	KK3=i/(N_CELL[0]*N_CELL[1]);
	x_L=KK1*Dx[0]+X_Low_E-(Delta_R_R)*sigma_max;
	x_H=x_L+Dx[0]+2*(Delta_R_R)*sigma_max;
	y_L=KK2*Dx[1]+Y_Low_E-(Delta_R_R)*sigma_max;
	y_H=y_L+Dx[1]+2*(Delta_R_R)*sigma_max;
	z_L=KK3*Dx[2]+Z_Low_E-(Delta_R_R)*sigma_max;
	z_H=z_L+Dx[2]+2*(Delta_R_R)*sigma_max;
	if ((((x_L>Inlet[i1].x1) && (x_L<Inlet[i1].x2)) | ((x_H>Inlet[i1].x1) && (x_H<Inlet[i1].x2)))|((x_L<Inlet[i1].x1) && (x_H>Inlet[i1].x2))) {
		if ((((y_L>Inlet[i1].y1) && (y_L<Inlet[i1].y2)) | ((y_H>Inlet[i1].y1) && (y_H<Inlet[i1].y2)))|((y_L<Inlet[i1].y1) && (y_H>Inlet[i1].y2))) {
			if ((((z_L>Inlet[i1].z1) && (z_L<Inlet[i1].z2)) | ((z_H>Inlet[i1].z1) && (z_H<Inlet[i1].z2)))|((z_L<Inlet[i1].z1) && (z_H>Inlet[i1].z2))) {
				return 1;
			} else {
				return 0;
			}
		} else {
			return 0;
		}
	} else {
		return 0;
	}
}
/***********************************************************************/
/***********************************************************************/
double Normal_Random(double m, double s)
/* ========================================================================
 * Returns a normal (Gaussian) distributed real number.
 * NOTE: use s > 0.0
 *
 * Uses a very accurate approximation of the normal idf due to Odeh & Evans, 
 * J. Applied Statistics, 1974, vol 23, pp 96-97.
 * ========================================================================
 */
{ 
  const double p0 = 0.322232431088;     const double q0 = 0.099348462606;
  const double p1 = 1.0;                const double q1 = 0.588581570495;
  const double p2 = 0.342242088547;     const double q2 = 0.531103462366;
  const double p3 = 0.204231210245e-1;  const double q3 = 0.103537752850;
  const double p4 = 0.453642210148e-4;  const double q4 = 0.385607006340e-2;
  double u, t, p, q, z;
  u   = Uniform_Random(0.000000001,.999999999999);
  if (u < 0.5)
    t = sqrt(-2.0 * log(u));
  else
    t = sqrt(-2.0 * log(1.0 - u));
  p   = p0 + t * (p1 + t * (p2 + t * (p3 + t * p4)));
  q   = q0 + t * (q1 + t * (q2 + t * (q3 + t * q4)));
  if (u < 0.5)
    z = (p / q) - t;
  else
    z = t - (p / q);
  return (m + s * z);
}
/***********************************************************************/
/***********************************************************************/
int Belong_Cell_to_Outlet(int i,int i1)  
{
	double x_L,y_L,z_L,x_H,y_H,z_H;
	int KK1,KK2,KK3;
	KK1=i%N_CELL[0];
	KK2=(i%(N_CELL[0]*N_CELL[1]))/N_CELL[0];
	KK3=i/(N_CELL[0]*N_CELL[1]);
	x_L=KK1*Dx[0]+X_Low_E-(Delta_R_R)*sigma_max;
	x_H=x_L+Dx[0]+2*(Delta_R_R)*sigma_max;
	y_L=KK2*Dx[1]+Y_Low_E-(Delta_R_R)*sigma_max;
	y_H=y_L+Dx[1]+2*(Delta_R_R)*sigma_max;
	z_L=KK3*Dx[2]+Z_Low_E-(Delta_R_R)*sigma_max;
	z_H=z_L+Dx[2]+2*(Delta_R_R)*sigma_max;
	R_Cut_Off=R_Cut_Off_R*sigma_max;
	if ((((x_L>Outlet[i1].x1) && (x_L<Outlet[i1].x2)) | ((x_H>Outlet[i1].x1) && (x_H<Outlet[i1].x2)))|((x_L<Outlet[i1].x1) && (x_H>Outlet[i1].x2))) {
		if ((((y_L>Outlet[i1].y1) && (y_L<Outlet[i1].y2)) | ((y_H>Outlet[i1].y1) && (y_H<Outlet[i1].y2)))|((y_L<Outlet[i1].y1) && (y_H>Outlet[i1].y2))) {
			if ((((z_L>Outlet[i1].z1) && (z_L<Outlet[i1].z2)) | ((z_H>Outlet[i1].z1) && (z_H<Outlet[i1].z2)))|((z_L<Outlet[i1].z1) && (z_H>Outlet[i1].z2))) {
				return 1;
			} else {
				return 0;
			}
		} else {
			return 0;
		}
	} else {
		return 0;
	}
}

