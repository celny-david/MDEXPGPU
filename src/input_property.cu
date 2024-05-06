/*
 ============================================================================
 Name        : input_property.cu
 Author      : klimam
 Version     :
 Description : Collects all parameters and functions relevant for initiation of simulation,
 	 	 	   seting its parameters and output properties

 NOTES : *
 ============================================================================
 */


#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "input_property.h"
#include "output.h"
#include "output_property.h"
#include "property.h"
#include "cuda_definition.h"
#include "lj_spline.h"
#include "lj_elspline.h"
#include "system_property.h"
#include "iterative_property.h"
#include "triplet.h"

/* Function: set_parameters_from_property
 * reads and set simulation properties from property.h
 *
 * Notes:
 * > The function modifies output_prop structure.
 * > Calls functions set_*
 *
 * BEWARE:
 * > CONSIDER THIS FUNCTIONALITY OBSOLETE -> FUTURE DELETION
 * Parameters:
 * void
 */
void set_parameters_from_property(void)
{
	set_molecule_count(p_mol_count,THREAD_P_BLOCK);
	set_mol_mass();
	set_l_from_density(p_density);

#if SUBSTANCE==LJF
	set_cutoff(3.0, 3.0*1.5);
#elif SUBSTANCE== ARGON
	set_cutoff(3.0*3.4, 3.0*3.4*1.5);
#elif SUBSTANCE== NITROGEN
	set_cutoff(9.9234, 9.9234*1.25);
#elif SUBSTANCE== SPCE
	set_cutoff(12.6622315608, 12.6875, 12.6875*1.5);
#else
	fprintf(stderr,"ERROR: set the substance cutoff and temperature");
	exit(-1);
#endif

	set_temperature(p_temperature);
	set_verlet_refresh(p_verlet_refresh);

	set_splines();
#if SUBSTANCE== SPCE
	set_elsplines();
#endif

	set_dt(p_d_dt);
	set_total_time(p_d_total_time);
	set_printout_steps(p_d_printout_steps);//for input -1 it does not print to file
	set_plb_hist_printout_steps(p_d_plb_hist_printout_steps);
	set_substance_name();
	set_output_name();
	set_n_shake(p_d_nshake);
	return;
}

/* Function: removeLastN
 * removes last n characters from given string and returns resulting string
 *
 * Parameters:
 * str - char* input string
 * n - int number of removed characters
 */
char* removeLastN(char* str, int n ) {
	int strLen = strlen(str);
    str[n <= strLen ? strLen-n : 0] = '\0';
    return str;
}

/* Function: set_input_from_file
 * sets input and output names and calls set_output_name
 *
 * Parameters:
 * input_file_whole - char* input string
 */
void set_input_from_file(char *input_file_whole)
{
	char name_buffer[256];

	inpp.input_source=1;
	strcpy(inpp.input_name,input_file_whole);
	if (strcmp(strrchr(inpp.input_name, '\0') - 5, ".cdef"))
	{
	printf("****ERROR**** \n");
	printf("Wrong input file suffix.\n");
	printf("Input file suffix must be .cdef !!\n");
	exit(-1);
	}
	strcpy(name_buffer,input_file_whole);
	printf("*** INFO - parsing parameters from : %s\n", name_buffer); // INFO so that the extension is also displayed
	set_output_name(removeLastN(name_buffer,5)); // NOTE wrong input handling is within the set_parameters_from_input_file function
	set_input_cfa_name(name_buffer);
	set_input_cfg_name(name_buffer);
	return;
}

/* Function: print_input_error
 * prints input error and exits program
 *
 * Parameters:
 * void
 */
void print_input_error(void)
{
	FILE * File_buffer;

	printf("***ERROR: wrong number of input parameters!\n");

	//File_buffer=fopen(outp.log_name,"a");
	//fprintf(File_buffer,"***ERROR: wrong number of input parameters!\n");
	//fclose(File_buffer);

	exit(-1);
}

/* Function: set_parameters_from_input_file
 * the function for setting the simulation properties from .def file
 *
 * Notes:
 * > calls set_outputs_and_periods function for output files generation settings
 *
 * Parameter:
 * no one, but reads from .def following parameters:
 *
 * N[0]= Number of requrested molecules (unsigned int) - this number is aligned to WORKLOAD
 * rho= Density of simulated system [kg*m^(-3)] (float)
 * T= System temperature [K] (float)
 * verlet_refresh= refreshing period of Verlet list [# steps] (unsigned int)
 * dt= integration step [ps] (float)
 * total_time= total simulation length [ps] (float)
 * n_shake= number of SHAKE iterations performed every steps [#]
 */
void set_parameters_from_input_file(void)
{
	FILE * Input_file;
	FILE * Printout_file;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
	// NOTE bit mask for what has been read 
	// unread means 0 bit at position
	
	uint read_mask=0; // initial setting as 0
	unsigned int ui_buffer=0;
	float fl_buffer=0.0;

	// NOTE temporary initialization until iterative schema is updated
	itep.plb_hist_printout_steps=1;
	itep.printout_steps=1;
	itep.printout=true;


	//Printout_file=fopen(output_prop.log_name,"w");
	//fclose(Printout_file);
	Printout_file=fopen(outp.log_name,"a");
	Input_file=fopen(inpp.input_name,"r");

	#if SUBSTANCE==LJF
		set_cutoff(3.0, 3.0*1.5);
	#elif SUBSTANCE== ARGON
		set_cutoff(3.0*3.4, 3.0*3.4*1.5);
	#elif SUBSTANCE== NITROGEN
		set_cutoff(9.9234, 9.9234*1.25);
	#elif SUBSTANCE== SPCE
		set_cutoff(12.6622315608, 12.6875, 12.6875*1.5);
	#else
		fprintf(stderr,"ERROR: set the substance cutoff and temperature");
		exit(-1);
	#endif

	set_mol_mass();
	set_substance_name();

	set_splines();
	// NOTE elsplines are empty when substance does not have lectrostatics enabled
	set_elsplines();
	
	if(Input_file==NULL)
	{
		printf("Error: Can not open the input file %s!\n",inpp.input_name);
		fprintf(Printout_file,"Error: Can not open the input file %s!\n",inpp.input_name);
		exit(-1);
	}

	while ((read = getline(&line, &len, Input_file)) != -1)
	{

		if(sscanf(line, "N[0]= %u\n",&ui_buffer)==1)
		{
			if (read_mask >>0 & 0x01)
			{
				printf("WARNING: duplicit input of N[0]. Latest taken.\n");
			}
			set_molecule_count(ui_buffer,THREAD_P_BLOCK);
			read_mask = read_mask | 1<<0;
		}

		if(sscanf(line, "rho= %f\n",&fl_buffer)==1)
		{
			if (read_mask >>1 & 0x01)
			{
				printf("WARNING: duplicit input of rho. Latest taken.\n");
			}
			set_l_from_density(fl_buffer);
			read_mask = read_mask | 1<<1;
		}

		if(sscanf(line, "T= %f\n",&fl_buffer)==1)
		{
			if (read_mask >>2 & 0x01)
			{
				printf("WARNING: duplicit input of T. Latest taken.\n");
			}
			set_temperature(fl_buffer);
			read_mask = read_mask | 1<<2;
		}

		if(sscanf(line, "verlet_refresh= %u\n",&ui_buffer)==1)
		{
			if (read_mask >>3 & 0x01)
			{
				printf("WARNING: duplicit input of verlet_refresh. Latest taken.\n");
			}
			set_verlet_refresh(ui_buffer);
			read_mask = read_mask | 1<<3;
		}


		if(sscanf(line, "dt= %f\n",&fl_buffer)==1)//integration step[ps]
		{
			if (read_mask >>4 & 0x01)
			{
				printf("WARNING: duplicit input of dt. Latest taken.\n");
			}
			set_dt(fl_buffer);
			read_mask = read_mask | 1<<4;
		}

		if(sscanf(line, "total_time= %f\n",&fl_buffer)==1)
		{
			if (read_mask >>5 & 0x01)
			{
				printf("WARNING: duplicit input of total_time. Latest taken.\n");
			}
			set_total_time(fl_buffer);
			read_mask = read_mask | 1<<5;
		}

		if(sscanf(line, "n_shake= %f\n",&fl_buffer)==1)
		{
			if (read_mask >>6 & 0x01)
			{
				printf("WARNING: duplicit input of n_shake. Latest taken.\n");
			}
			set_n_shake(fl_buffer);
			read_mask = read_mask | 1<<6;
		}
		if(sscanf(line, "l_x= %f\n",&fl_buffer)==1)
		{
			if (read_mask >>7 & 0x01)
			{
				printf("WARNING: duplicit input of l_x. Latest taken.\n");
			}
			sysp.lx = (value_t) fl_buffer;
			read_mask = read_mask | 1<<7;
		}
		if(sscanf(line, "l_y= %f\n",&fl_buffer)==1)
		{
			if (read_mask >>8 & 0x01)
			{
				printf("WARNING: duplicit input of l_y. Latest taken.\n");
			}
			sysp.ly = (value_t) fl_buffer;
			read_mask = read_mask | 1<<8;
		}
		if(sscanf(line, "l_z= %f\n",&fl_buffer)==1)
		{
			if (read_mask >>9 & 0x01)
			{
				printf("WARNING: duplicit input of l_z. Latest taken.\n");
			}
			sysp.lz = (value_t) fl_buffer;
			read_mask = read_mask | 1<<9;
		}		
	}
	free(line);

	// Ready for update with unreaded variable specification

	// DEBUG show the bit mask
	// int bit;
	// for(bit = (sizeof(unsigned int) * 8)-1; bit>=0; bit--)
	// {
	// 	printf("%2i ",bit);        
	// }
	// printf("\n-----------------------------------------------------------------------------------------------\n");
	// for(bit = (sizeof(uint) * 8)-1; bit>=0; bit--)
	// {
	// 	printf("%2i ", (read_mask & 1<<bit) >>bit & 0x01);
	// }
	// printf("\n");
	//DEBUG END

	fclose(Input_file);
	set_outputs_and_periods();

	if ( !(read_mask >>0 & 0x01))
	{
		printf("Can't read number of molecules N[0]!!\n");
		fprintf(Printout_file,"Can't read number of molecules N[0]!!\n");
		exit(-1);
	}
	if ( !(read_mask >>1 & 0x01) && !(read_mask >>7 & 0x0111))
	{
		printf("Can't read density rho and not all lx,ly,lz are given!!\n");
		fprintf(Printout_file,"Can't read density rho and not all lx,ly,lz are given!!\n");
		exit(-1);
	}
	if ( !(read_mask >>2 & 0x01))
	{
		printf("Can't read temperature T!!\n");
		fprintf(Printout_file,"Can't read temperature T!!\n");
		exit(-1);
	}
	if ( !(read_mask >>3 & 0x01))
	{
		printf("Can't read verlet neighbour list size verlet_neighbour!!\n");
		fprintf(Printout_file,"Can't read verlet neighbour list size verlet_neighbour!!\n");
		exit(-1);
	}
	if ( !(read_mask >>4 & 0x01))
	{
		printf("Can't read integration step length dt!!\n");
		fprintf(Printout_file,"Can't read integration step length dt!!\n");
		exit(-1);
	}
	if ( !(read_mask >>5 & 0x01))
	{
		printf("Can't read total_time!!\n");
		fprintf(Printout_file,"Can't read total_time!!\n");
		exit(-1);
	}
	if ( !(read_mask >>6 & 0x01))
	{
		printf("Can't read n_shake!!\n");
		fprintf(Printout_file,"Can't read n_shake!!\n");
		exit(-1);
	}
	if ( !(read_mask >>1 & 0x01) && !(read_mask >>7 & 0x01))
	{
		printf("Can't read l_x, and rho is not given!!\n");
		fprintf(Printout_file,"Can't read l_x, and rho is not given!!\n");
		exit(-1);
	}
	if ( !(read_mask >>1 & 0x01) && !(read_mask >>8 & 0x01))
	{
		printf("Can't read l_y, and rho is not given!!\n");
		fprintf(Printout_file,"Can't read l_y, and rho is not given!!\n");
		exit(-1);
	}
	if ( !(read_mask >>1 & 0x01) && !(read_mask >>9 & 0x01))
	{
		printf("Can't read l_z, and rho is not given!!\n");
		fprintf(Printout_file,"Can't read l_z, and rho is not given!!\n");
		exit(-1);
	}
	fclose(Printout_file);

// TODO check the sense of this here
#ifdef EXPANSION
	set_box_file_name();
	L_from_density_file(0);
#endif
	return;
}

/* Function: set_outputs_and_periods
 * the function for setting output files and its writing period from .def file
 *
 * Parameter:
 * no one, but reads from .def following parameters:
 *
 * dt.all=period[ps], write all possible files with printout period [ps]
 * dt.mac=period[ps], write MACSIMUS standart files (.plb and .cp) with printout period [ps]
 *
 * dt.plb=period[ps], write .plb file with printout period [ps]
 * same is for:
 * dt.cp, dt.cfg, dt.cfa, dt.test, dt.test2, dt.test3, dt.hist, dt.eng
 *
 */
void set_outputs_and_periods(void)
{
	float fl_buffer=0;
	unsigned int  ui_buffer=0;
	FILE *Input_file;
	char * line = NULL;
	size_t len = 0;
    ssize_t read;
    Input_file=fopen(inpp.input_name,"r");

	while ((read = getline(&line, &len, Input_file)) != -1) {

		if(sscanf(line, "dt.plb= %f\n",&fl_buffer)==1)
		{
			ui_buffer=ceil(fl_buffer/itep.dt);
			outp.plb_out=true;
			outp.plb_print_period=ui_buffer;
		}

		if(sscanf(line, "dt.cp= %f\n",&fl_buffer)==1)
		{
			ui_buffer=ceil(fl_buffer/itep.dt);
			outp.cp_out=true;
			outp.cp_print_period=ui_buffer;
		}

		if(sscanf(line, "dt.cfg= %f\n",&fl_buffer)==1)
		{
			ui_buffer=ceil(fl_buffer/itep.dt);
			outp.cfg_out=true;
			outp.cfg_print_period=ui_buffer;
		}

		if(sscanf(line, "dt.cfa= %f\n",&fl_buffer)==1)
		{
			ui_buffer=ceil(fl_buffer/itep.dt);
			outp.cfa_out=true;
			outp.cfa_print_period=ui_buffer;
		}
		
		if(sscanf(line, "dt.vel= %f\n",&fl_buffer)==1)
		{
			ui_buffer=ceil(fl_buffer/itep.dt);
			outp.vel_out=true;
			outp.vel_print_period=ui_buffer;
		}

		if(sscanf(line, "dt.oob= %f\n",&fl_buffer)==1)
		{
			ui_buffer=ceil(fl_buffer/itep.dt);
			outp.oob_out=true;
			outp.oob_print_period=ui_buffer;
		}

		if(sscanf(line, "dt.test4= %f\n",&fl_buffer)==1)
		{
			ui_buffer=ceil(fl_buffer/itep.dt);
			outp.test4_out=true;
			outp.test4_print_period=ui_buffer;
		}

		if(sscanf(line, "dt.fcs= %f\n",&fl_buffer)==1)
		{
			ui_buffer=ceil(fl_buffer/itep.dt);
			outp.fcs_out=true;
			outp.fcs_print_period=ui_buffer;
		}

		if(sscanf(line, "dt.hist= %f\n",&fl_buffer)==1)
		{
			ui_buffer=ceil(fl_buffer/itep.dt);
			outp.history_out=true;
			outp.history_print_period=ui_buffer;
		}

		if(sscanf(line, "dt.eng= %f\n",&fl_buffer)==1)
		{
			ui_buffer=ceil(fl_buffer/itep.dt);
			outp.energy_out=true;
			outp.energy_print_period=ui_buffer;
		}

		if(sscanf(line, "dt.all= %f\n",&fl_buffer)==1)
		{
			ui_buffer=ceil(fl_buffer/itep.dt);
			outp.plb_out=true;
			outp.cp_out=true;
			outp.cfg_out=true;
			outp.cfa_out=true;
			outp.fcs_out=true;
			outp.vel_out=true;
			outp.oob_out=true;
			outp.test4_out=true;
			outp.history_out=true;
			outp.energy_out=true;
			outp.plb_print_period=ui_buffer;
			outp.cp_print_period=ui_buffer;
			outp.cfg_print_period=ui_buffer;
			outp.cfa_print_period=ui_buffer;
			outp.fcs_print_period=ui_buffer;
			outp.vel_print_period=ui_buffer;
			outp.oob_print_period=ui_buffer;
			outp.test4_print_period=ui_buffer;
			outp.history_print_period=ui_buffer;
			outp.energy_print_period=ui_buffer;
			set_printout_steps(ui_buffer);//Remove after iterative schema update
		}

		if(sscanf(line, "dt.mac= %f\n",&fl_buffer)==1)
		{
			ui_buffer=ceil(fl_buffer/itep.dt);
			outp.plb_out=true;
			outp.cp_out=true;
			outp.plb_print_period=ui_buffer;
			outp.cp_print_period=ui_buffer;
			set_printout_steps(ui_buffer); //Remove after iterative schema update

		}
	}
	fclose(Input_file);
	free(line);
	return;
}

/* Function: L_from_density_file
 * reads density-time dependency from simname.box[kg/m3] file, linearly interpolates
 * an returns size[Angstrom]
 *
 * Parameter:
 * t = value_t simulation time[picosecond]
 *
 */
void L_from_density_file(value_t t) /******************************** Lfromfile */
{
  typedef double vector[3];

  static struct history_s {
    double t;
    vector L;
  } one,hist[3];
  int k=0;
  static int nl=0;
  static FILE *rhof=NULL;
  char line[128];
  double rho,boxL;
  int i;
  //char *box_file;

  //box_file=removeLastN(inpp.input_name,4);
  //strcat (box_file,"box");
  //printf("%s\n",box_file);
 // if (!rhof) {
    rhof=fopen(inpp.box_file_name,"rt");
   // }
while (t>=hist[2].t || nl<3) {
    if (!fgets(line,128,rhof)) break;
    if (-1==-1) {
      if (sscanf(line,"%lf%lf",&one.t,&rho)<2) continue;
      one.L[0]=one.L[1]=one.L[2]=cbrt(sysp.system_real_mass /rho);}
    else
      if (sscanf(line,"%lf%lf%lf%lf",&one.t,one.L,one.L+1,one.L+2)<4) continue;
    hist[0]=hist[1]; hist[1]=hist[2]; hist[2]=one;
    nl++; }

  i=t>hist[1].t;
  fclose(rhof);
  boxL= 0.1e11*(hist[i].L[k]+(t-hist[i].t)/(hist[i+1].t-hist[i].t)*(hist[i+1].L[k]-hist[i].L[k]));
  set_l(boxL,boxL,boxL);
  set_density_from_l(boxL);
  return;
}

void set_box_file_name(void)
{
	char *box_file;
	box_file=removeLastN(inpp.input_name,4);
	strcat(box_file,"box");
	strcpy(inpp.box_file_name,box_file);
	return;
}
/* Function: set_input_cfa_name
 * reads and set name of input .cfa file containing initial configuration
 *
 * Parameters:
 * input_cfa_name_buffer -char[] name of input .cfa file
 */
void set_input_cfa_name(char input_cfa_name_buffer[])
{

	if(strlen(input_cfa_name_buffer)>=(256-4) )
	{
		fprintf (stderr,"Error input_cfa name length exceeded allowed limit!\n");
		return;
	}
	else
	{
	strcpy(inpp.input_cfa_name,input_cfa_name_buffer);
	strcat(inpp.input_cfa_name,".cfa");
	}

	return;
}

/* Function: set_input_cfg_name
 * reads and set name of input .cfg file containing initial configuration
 *
 * Parameters:
 * input_cfg_name_buffer -char[] name of input .cfg file
 */
void set_input_cfg_name(char input_cfg_name_buffer[])
{

	if(strlen(input_cfg_name_buffer)>=(256-4) )
	{
		fprintf (stderr,"Error input_cfg name length exceeded allowed limit!\n");
		return;
	}
	else
	{
	strcpy(inpp.input_cfg_name,input_cfg_name_buffer);
	strcat(inpp.input_cfg_name,".cfg");
	}

	return;
}

/* Function: set_parameters_from_cdef_cfa_file
 * the function for setting the simulation properties from .def file and parses initial configuration(positions and velocities) from .cfa file
 *
 * Notes:
 * > calls set_outputs_and_periods function for output files generation settings
 * > calls position_decombiner to obtain relative atomar and molecular centers of gravity coordinates and velocities
 *
 * Parameter:
 * molecule_position *triplet coordinates of molecular centers of gravity
 * atom_position *triplet relative atom coordinates
 * molecule_velocity *triplet velocities of molecular centers of gravity
 * atom_velocity *triplet relative atom velocities
 * atom_force *triplet forces acting on atom, set on 0 by default
 *
 *
 * N[0]= Number of molecules in simulation (unsigned int)
 * rho= Density of simulated system [kg*m^(-3)] (float)
 * T= System temperature [K] (float)
 * verlet_neighbour= length of Verlet list (unsigned int)
 * verlet_refresh= refreshing period of Verlet list [steps???] (unsigned int)
 * dt= integration step [ps] (float)
 * total_time= total simulation length [ps] (float)
 * n_shake= number of SHAKE iterations
 *
 *
 */
void set_parameters_from_cdef_cfa_file(triplet *molecule_position, triplet *atom_position,
									   triplet *molecule_velocity, triplet *atom_velocity,
	    							   triplet *atom_force)
{
	FILE * Input_file;
	FILE * Printout_file;
	FILE * CFA_file;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;

	//itep.plb_hist_printout_steps=1;
	//itep.printout_steps=1;
	//itep.printout=true;

	triplet *input_velocity=(triplet*) malloc(sizeof(triplet));
	triplet *input_position=(triplet*) malloc(sizeof(triplet));

	int read_mark=0; // NOTE better check is with the bit masking -> allow detection of what failed
	unsigned int ui_buffer=0,i=0,ui_buffer2=0;
	double fl_buffer=0.0,fl_buffer2=0.0;
	double buffer[6];
	unsigned int Inputs_number=4;

	//Printout_file=fopen(output_prop.log_name,"w");
	//fclose(Printout_file);
	Printout_file=fopen(outp.log_name,"a");
	Input_file=fopen(inpp.input_name,"r");
	CFA_file=fopen(inpp.input_cfa_name,"r");

	#if SUBSTANCE==ARGON
		set_cutoff(3.0, 3.0*1.5);
	#elif SUBSTANCE== NITROGEN
		set_cutoff(9.9234, 9.9234*1.25);
	#elif SUBSTANCE== SPCE
		set_cutoff(12.6622315608, 12.6875, 12.6875*1.5);
	#else
		fprintf(stderr,"ERROR: set the substance cutoff and temperature");
		exit(-1);
	#endif

	set_mol_mass();
	set_substance_name();

	set_splines();
	#if SUBSTANCE== SPCE
		set_elsplines();
	#endif

	if(Input_file==NULL)
	{
		printf("Error: Can not open the input file %s!\n",inpp.input_name);
		fprintf(Printout_file,"Error: Can not open the input file %s!\n",inpp.input_name);
		exit(-1);
	}
	if(CFA_file==NULL)
	{
		printf("Error: Can not open the cfa input file %s!\n",inpp.input_cfa_name);
		fprintf(Printout_file,"Error: Can not open the cfa input file %s!\n",inpp.input_cfa_name);
		exit(-1);
	}


	while ((read = getline(&line, &len, CFA_file)) != -1)
	{

		if(strstr(line,"N ns (of species 0)")!=NULL)
		{
		  sscanf(line, "%u %u",&ui_buffer,&ui_buffer2);
		  sysp.molecule_count = ui_buffer;
		  sysp.atom_per_molecule =SUBSTANCE;
		  sysp.atom_count = sysp.molecule_count*sysp.atom_per_molecule;
		  if(ui_buffer2 != SUBSTANCE)
		  {
			  printf("ERROR .cfa file was generated for another substance than mac_module compiled for!\n ");
			  printf("mac_module atom_per_molecule: %u!\n ",sysp.atom_per_molecule);
			  printf(".cfa file atom_per_molecule: %u!\n ",ui_buffer2);
			  fprintf(Printout_file,"ERROR .cfa file was generated for another substance than mac_module compiled for!\n");
			  fprintf(Printout_file,"mac_module atom_per_molecule: %u!\n ",sysp.atom_per_molecule);
			  fprintf(Printout_file,".cfa file atom_per_molecule: %u!\n ",ui_buffer2);
			  exit(-1);

		  }
		  else if(sysp.atom_count%THREAD_P_BLOCK==0)
		  {
			  sysp.molecule_count_aligned = sysp.molecule_count;
			  sysp.atom_count_aligned = sysp.molecule_count*sysp.atom_per_molecule;
		  }
		  else
		  {
			  printf("ERROR unaligned molecule count in .cdef file!\n");
			  fprintf(Printout_file,"ERROR unaligned molecule count in .cdef file!\n\n");
			  exit(-1);
		  }
		}



		if(strstr(line,"key ns L[3]")!=NULL)
		{

		  sscanf(line, "%u %u %lf %lf %lf",&ui_buffer,&ui_buffer,&buffer[0],&buffer[1],&buffer[2]);
		  sysp.lx=buffer[0];
		  sysp.ly=buffer[1];
		  sysp.lz=buffer[2];

		  set_density();
		}


		if(strstr(line,"En.tot")!=NULL)
		{

		  sscanf(line, "%u %u %lf %lf",&ui_buffer,&ui_buffer,&fl_buffer2,&fl_buffer);
		  itep.dt=fl_buffer;
		}



		if(i>40)//reads just first 40 lines
		{
			break;
		}
		i++;
	}

	printf("*** INFO - parsing configuration from : %s\n", inpp.input_cfa_name);
	fprintf(Printout_file,"*** INFO - parsing configuration from : %s\n",inpp.input_cfa_name);

	// = allocation section =
	triplet_alloc(molecule_position,sysp.molecule_count_aligned);
	triplet_alloc(atom_position,sysp.atom_count_aligned);

	triplet_alloc(molecule_velocity,sysp.molecule_count_aligned);
	triplet_alloc(atom_velocity,sysp.atom_count_aligned);

	triplet_alloc(atom_force,sysp.atom_count_aligned);

	triplet_alloc(input_position,sysp.atom_count_aligned);
	triplet_alloc(input_velocity,sysp.atom_count_aligned);

	i=0;
	while ((read = getline(&line, &len, CFA_file)) != -1)
	{
		if(strstr(line,"site[")!=NULL)
		{
		  sscanf(line, " %lf %lf %lf %lf %lf %lf",&buffer[0],&buffer[1],&buffer[2],&buffer[3],&buffer[4],&buffer[5]);
		  input_position->x[i]=buffer[0];
		  input_position->y[i]=buffer[1];
		  input_position->z[i]=buffer[2];
		  input_velocity->x[i]=buffer[3];
		  input_velocity->y[i]=buffer[4];
		  input_velocity->z[i]=buffer[5];
		  i++;
		}
	}

	fclose(CFA_file);
	free(line);

	while ((read = getline(&line, &len, Input_file)) != -1)
	{

		if(sscanf(line, "verlet_refresh= %u\n",&ui_buffer)==1)
		{
			set_verlet_refresh(ui_buffer);
			read_mark++;
		}

		if(sscanf(line, "total_time= %lf\n",&fl_buffer)==1)
		{
			set_total_time(fl_buffer);
			read_mark++;
		}

		if(sscanf(line, "n_shake= %lf\n",&fl_buffer)==1)
		{
			set_n_shake(fl_buffer);
			read_mark++;
		}
	}
	if(read_mark!=Inputs_number)
	{
		printf("Can't read some input variable from .cdef file, .cfa file exists!\n");
		fprintf(Printout_file,"Can't read some input variable from .cdef file, .cfa file exists!\n");
	}

	if(read_mark!=Inputs_number)
	{
		exit(-1);
	}

	fclose(Input_file);
	set_outputs_and_periods();
    set_max_step();
	position_decombiner(input_position,molecule_position,atom_position);
	position_decombiner(input_velocity,molecule_velocity,atom_velocity);

	for(unsigned int i=0;i<sysp.atom_count_aligned;i++)
	{
		atom_force->x[i] = 0.0;
		atom_force->y[i] = 0.0;
		atom_force->z[i] = 0.0;
	}
	set_temperature(Return_Temperature(input_velocity));
	triplet_dealloc(input_position);
	triplet_dealloc(input_velocity);

	fclose(Printout_file);

	return;
}

/* Function: set_parameters_from_cdef_cfg_file
 * the function for setting the simulation properties from .def file and parses initial configuration(positions and velocities) from .cfg binary file
 *
 * Notes:
 * > calls set_outputs_and_periods function for output files generation settings
 * > calls position_decombiner to obtain relative atomar and molecular centers of gravity coordinates and velocities
 *
 * Parameter:
 * molecule_position *triplet coordinates of molecular centers of gravity
 * atom_position *triplet relative atom coordinates
 * molecule_velocity *triplet velocities of molecular centers of gravity
 * atom_velocity *triplet relative atom velocities
 * atom_force *triplet forces acting on atom, set on 0 by default
 *
 *
 * N[0]= Number of molecules in simulation (unsigned int)
 * rho= Density of simulated system [kg*m^(-3)] (float)
 * T= System temperature [K] (float)
 * verlet_neighbour= length of Verlet list (unsigned int)
 * verlet_refresh= refreshing period of Verlet list [steps???] (unsigned int)
 * dt= integration step [ps] (float)
 * total_time= total simulation length [ps] (float)
 * n_shake= number of SHAKE iterations
 *
 */
void set_parameters_from_cdef_cfg_file(triplet *molecule_position, triplet *atom_position,
									   triplet *molecule_velocity, triplet *atom_velocity,
	    							   triplet *atom_force)
{
	FILE * Input_file;
	FILE * Printout_file;
	FILE * CFG_file;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;

	triplet *input_velocity=(triplet*) malloc(sizeof(triplet));
	triplet *input_position=(triplet*) malloc(sizeof(triplet));

	int read_mark=0; // NOTE better check is with the bit masking -> allow detection of what failed
	unsigned int ui_buffer=0,ui_buffer2=0;
	double fl_buffer=0.0;
	unsigned int Inputs_number=4;
    int bufrInt3[3], bufrOptions[34], bufrInt4[4], bufrInt,bufrInt2[2];
    double bufrDoub3[3],bufrDoub;



	//Printout_file=fopen(output_prop.log_name,"w");
	Printout_file=fopen(outp.log_name,"a");
	Input_file=fopen(inpp.input_name,"r");
	CFG_file=fopen(inpp.input_cfg_name,"rb");

	#if SUBSTANCE==ARGON
		set_cutoff(3.0, 3.0*1.5);
	#elif SUBSTANCE== NITROGEN
		set_cutoff(9.9234, 9.9234*1.25);
	#elif SUBSTANCE== SPCE
		set_cutoff(12.6622315608, 12.6875, 12.6875*1.5);
	#else
		fprintf(stderr,"ERROR: set the substance cutoff and temperature");
		exit(-1);
	#endif

	set_mol_mass();
	set_substance_name();

	set_splines();
	#if SUBSTANCE== SPCE
		set_elsplines();
	#endif

	if(Input_file==NULL)
	{
		printf("Error: Can not open the input file %s!\n",inpp.input_name);
		fprintf(Printout_file,"Error: Can not open the input file %s!\n",inpp.input_name);
		exit(-1);
	}
	if(CFG_file==NULL)
	{
		printf("Error: Can not open the cfg input file %s!\n",inpp.input_cfg_name);
		fprintf(Printout_file,"Error: Can not open the cfg input file %s!\n",inpp.input_cfg_name);
		exit(-1);
	}

	  fread(&bufrInt3,sizeof(bufrInt3),1,CFG_file);
	  //printf("cfgkey: %d \n",bufrInt3[1]);
	  fread(&bufrOptions,sizeof(bufrOptions),1,CFG_file);
	  //printf("options[32]: %d \n",bufrOptions[32]);
	  fread(&bufrInt3,sizeof(bufrInt3),1,CFG_file);
	  //printf("nspec: %d \n",bufrInt3[1]);

	  fread(&bufrInt4,sizeof(bufrInt4),1,CFG_file);
	  //printf("N ns (of species 0) %d %d\n",bufrInt4[1],bufrInt4[2]);

	  sysp.molecule_count = bufrInt4[1];
	  sysp.atom_per_molecule =SUBSTANCE;
	  sysp.atom_count = sysp.molecule_count*sysp.atom_per_molecule;
	  if(bufrInt4[2] != SUBSTANCE)
	  {
		  printf("ERROR .cfg file was generated for another substance than mac_module compiled for!\n ");
		  printf("mac_module atom_per_molecule: %u!\n ",sysp.atom_per_molecule);
		  printf(".cfg file atom_per_molecule: %u!\n ",ui_buffer2);
		  fprintf(Printout_file,"ERROR .cfg file was generated for another substance than mac_module compiled for!\n");
		  fprintf(Printout_file,"mac_module atom_per_molecule: %u!\n ",sysp.atom_per_molecule);
		  fprintf(Printout_file,".cfg file atom_per_molecule: %u!\n ",ui_buffer2);
		  exit(-1);

	  }
	  else if(sysp.atom_count%THREAD_P_BLOCK==0)
	  {
		  sysp.molecule_count_aligned = sysp.molecule_count;
		  sysp.atom_count_aligned = sysp.molecule_count*sysp.atom_per_molecule;
	  }
	  else
	  {
		  printf("ERROR unaligned molecule count in .cdef file!\n");
		  fprintf(Printout_file,"ERROR unaligned molecule count in .cdef file!\n\n");
		  exit(-1);
	  }

	  fread(&bufrInt3,sizeof(bufrInt3),1,CFG_file);
      fread(&bufrDoub3,sizeof(bufrDoub3),1,CFG_file);
	  fread(&bufrInt,sizeof(bufrInt),1,CFG_file);
	  //printf("key N t h En.tot %d %d %lf %lf %lf\n",bufrInt3[1],bufrInt3[2],bufrDoub3[0],bufrDoub3[1],bufrDoub3[2]);
	  itep.dt=bufrDoub3[1];


	  fread(&bufrInt3,sizeof(bufrInt3),1,CFG_file);
	  fread(&bufrDoub3,sizeof(bufrDoub3),1,CFG_file);
	  fread(&bufrInt,sizeof(bufrInt),1,CFG_file);
	  //printf("key ns L[3] %d %d %lf %lf %lf\n",bufrInt3[1],bufrInt3[2],bufrDoub3[0],bufrDoub3[1],bufrDoub3[2]);
	  sysp.lx=bufrDoub3[0];
	  sysp.ly=bufrDoub3[1];
	  sysp.lz=bufrDoub3[2];
	  set_density();


	printf("*** INFO - parsing configuration from : %s\n", inpp.input_cfg_name);
	fprintf(Printout_file,"*** INFO - parsing configuration from : %s\n",inpp.input_cfg_name);

	// = allocation section =
	triplet_alloc(molecule_position,sysp.molecule_count_aligned);
	triplet_alloc(atom_position,sysp.atom_count_aligned);

	triplet_alloc(molecule_velocity,sysp.molecule_count_aligned);
	triplet_alloc(atom_velocity,sysp.atom_count_aligned);

	triplet_alloc(atom_force,sysp.atom_count_aligned);

	triplet_alloc(input_position,sysp.atom_count_aligned);
	triplet_alloc(input_velocity,sysp.atom_count_aligned);



	fread(&bufrInt3,sizeof(bufrInt3),1,CFG_file);
	fread(&bufrDoub3,sizeof(bufrDoub3),1,CFG_file);
	fread(&bufrInt,sizeof(bufrInt),1,CFG_file);
	//printf("key intval vecval[3] %d %d %lf %lf %lf\n",bufrInt3[1],bufrInt3[2],bufrDoub3[0],bufrDoub3[1],bufrDoub3[2]);

	fread(&bufrInt,sizeof(bufrInt),1,CFG_file);
	fread(&bufrDoub,sizeof(bufrDoub),1,CFG_file);
	//printf("unknown %lf\n",bufrDoub);
	fread(&bufrDoub,sizeof(bufrDoub),1,CFG_file);
	//printf("unknown %lf\n",bufrDoub);
	fread(&bufrDoub3,sizeof(bufrDoub3),1,CFG_file);
	//printf("a[0]lambda[0..3] %lf %lf %lf\n",bufrDoub3[0],bufrDoub3[1],bufrDoub3[2]);

	fread(&bufrDoub3,sizeof(bufrDoub3),1,CFG_file);
	//printf("unknown: %lf %lf %lf\n",bufrDoub3[0],bufrDoub3[1],bufrDoub3[2]);

	for(unsigned int i=0;i<sysp.atom_count_aligned;i++)
	{
	  fread(&bufrDoub3,sizeof(bufrDoub3),1,CFG_file);
	  //  printf("xyz position of %u th atom (site[%u]) %lf %lf %lf\n",i,i,bufrDoub3[0],bufrDoub3[1],bufrDoub3[2]);
	  input_position->x[i]=bufrDoub3[0];
	  input_position->y[i]=bufrDoub3[1];
	  input_position->z[i]=bufrDoub3[2];
	}

	fread(&bufrInt2,sizeof(bufrInt2),1,CFG_file);
	fread(&bufrDoub,sizeof(bufrDoub),1,CFG_file);
	//printf("unknown %lf\n",bufrDoub);
	fread(&bufrDoub,sizeof(bufrDoub),1,CFG_file);
	//printf("unknown %lf\n",bufrDoub);
	fread(&bufrDoub3,sizeof(bufrDoub3),1,CFG_file);
	//printf("unknown %lf %lf %lf\n",bufrDoub3[0],bufrDoub3[1],bufrDoub3[2]);
	fread(&bufrDoub3,sizeof(bufrDoub3),1,CFG_file);
	//printf("unknown %lf %lf %lf\n",bufrDoub3[0],bufrDoub3[1],bufrDoub3[2]);

	for(unsigned int i=0;i<sysp.atom_count_aligned;i++)
	{
	  fread(&bufrDoub3,sizeof(bufrDoub3),1,CFG_file);
	  //printf("xyz velocity of %u th atom (site[%u]) %lf %lf %lf\n",i,i,bufrDoub3[0],bufrDoub3[1],bufrDoub3[2]);
	  input_velocity->x[i]=bufrDoub3[0];
	  input_velocity->y[i]=bufrDoub3[1];
	  input_velocity->z[i]=bufrDoub3[2];
	}


	fclose(CFG_file);
	free(line);

	while ((read = getline(&line, &len, Input_file)) != -1)
	{
		
		if(sscanf(line, "verlet_refresh= %u\n",&ui_buffer)==1)
		{
			set_verlet_refresh(ui_buffer);
			read_mark++;
		}

		if(sscanf(line, "total_time= %lf\n",&fl_buffer)==1)
		{
			set_total_time(fl_buffer);
			read_mark++;
		}

		if(sscanf(line, "n_shake= %lf\n",&fl_buffer)==1)
		{
			set_n_shake(fl_buffer);
			read_mark++;
		}
	}
	if(read_mark!=Inputs_number)
	{
		printf("Can't read some input variable from .cdef file, .cfg file exists!\n");
		fprintf(Printout_file,"Can't read some input variable from .cdef file, .cfg file exists!\n");
	}

	if(read_mark!=Inputs_number)
	{
		exit(-1);
	}

	fclose(Input_file);
	set_outputs_and_periods();
	set_max_step();
	position_decombiner(input_position,molecule_position,atom_position);
	position_decombiner(input_velocity,molecule_velocity,atom_velocity);

	for(unsigned int i=0;i<sysp.atom_count_aligned;i++)
	{
		atom_force->x[i] = 0.0;
		atom_force->y[i] = 0.0;
		atom_force->z[i] = 0.0;
	}
	set_temperature(Return_Temperature(input_velocity));
	triplet_dealloc(input_position);
	triplet_dealloc(input_velocity);

	fclose(Printout_file);

	return;
}

/* Function: position_decombiner
 * splits global coordinates into coordinates of molecular center of gravity and relative atom coordinates
 *
 * Notes:
 * > can be used for velocities too
 *
 * Parameter:
 * position_in *triplet input global coordinates
 * molecule_postition_in *triplet output recalculated coordinates of molecular centers of gravity
 * atom_position_in *triplet output recalculated relative atom coordinates
 *
 *
 */
void position_decombiner(triplet *position_in,triplet *molecule_postition_in,triplet *atom_position_in)
{

	for(unsigned int i=0;i<sysp.molecule_count_aligned;i++)
	{
		molecule_postition_in->x[i]=0.0f;
		molecule_postition_in->y[i]=0.0f;
		molecule_postition_in->z[i]=0.0f;
		for(unsigned int k=0;k<sysp.atom_per_molecule;k++)
		{
			molecule_postition_in->x[i]+=sysp.mol_mass_template[k]*position_in->x[sysp.atom_per_molecule*i+k];
			molecule_postition_in->y[i]+=sysp.mol_mass_template[k]*position_in->y[sysp.atom_per_molecule*i+k];
			molecule_postition_in->z[i]+=sysp.mol_mass_template[k]*position_in->z[sysp.atom_per_molecule*i+k];
		}
		molecule_postition_in->x[i]*=sysp.mol_mass_template[sysp.atom_per_molecule];
		molecule_postition_in->y[i]*=sysp.mol_mass_template[sysp.atom_per_molecule];
		molecule_postition_in->z[i]*=sysp.mol_mass_template[sysp.atom_per_molecule];

		for(unsigned int k=0;k<sysp.atom_per_molecule;k++)
		{
			atom_position_in->x[sysp.atom_per_molecule*i+k]=position_in->x[sysp.atom_per_molecule*i+k]-molecule_postition_in->x[i];
			atom_position_in->y[sysp.atom_per_molecule*i+k]=position_in->y[sysp.atom_per_molecule*i+k]-molecule_postition_in->y[i];
			atom_position_in->z[sysp.atom_per_molecule*i+k]=position_in->z[sysp.atom_per_molecule*i+k]-molecule_postition_in->z[i];
		}
	}

}

value_t get_lx_from_density_file(value_t t, FILE* rhof)
{
  typedef double vector[3];

  static struct history_s {
    double t;
    vector L;
  } one,hist[3];

  int k=0;
  static int nl=0;
  char line[128];
  double rho;
  int i;

  while (t>=hist[2].t || nl<3)
  {
    if (!fgets(line,128,rhof)) break;
    if (-1==-1) {
      if (sscanf(line,"%lf%lf",&one.t,&rho)<2) continue;
      one.L[0]=one.L[1]=one.L[2]=cbrt(sysp.system_real_mass /rho);}
    else
      if (sscanf(line,"%lf%lf%lf%lf",&one.t,one.L,one.L+1,one.L+2)<4) continue;
    hist[0]=hist[1]; hist[1]=hist[2]; hist[2]=one;
    nl++;
  }

  i=t>hist[1].t;
  return 0.1e11*(hist[i].L[k]+(t-hist[i].t)/(hist[i+1].t-hist[i].t)*(hist[i+1].L[k]-hist[i].L[k]));
}

void fill_lx_buffer(value_t *lx_buffer,unsigned int buffer_size)
{
	FILE * BoxFile;
	BoxFile=fopen(inpp.box_file_name,"rt");

	for(unsigned int i=0;i<buffer_size;i++)
	{
		lx_buffer[i]=get_lx_from_density_file(i*0.5*itep.dt, BoxFile);
	}

	fclose(BoxFile);
}
