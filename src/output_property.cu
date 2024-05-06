/*
 * output_property.cu
 *
 *  Created on: Jan 23, 2020
 *      Author: dixiw
 */

#include <stdio.h>
#include <string.h> // strcat

#include "output_property.h"
#include "system_property.h"
#include "iterative_property.h"
#include "constants.h"


void set_output_name(char output_name_buffer[])
{

	if(strlen(output_name_buffer)>=(256-4) )
	{
		fprintf (stderr,"Error output_name_buffer length exceeded allowed limit!\n");
		return;
	}
	else
	{
	strcpy(outp.output_files_name,output_name_buffer);
	}

	return;
}

void set_output_name(void)
{
	sprintf(outp.output_files_name,"%s_alg%i_type%i_%s_t%i_%04i",sysp.name,3,POT_TYPE,(PREC == FLOAT) ? "fl" : "db",(int)(1+itep.total_time/(1.0*itep.dt)),int(itep.dt*1000));
	return;
}

void create_files(void)
{
	FILE* File_buffer;

	if (outp.fcs_out)
	{
		strcpy(outp.fcs_name,outp.output_files_name);
		strcat(outp.fcs_name,".fcs");
		File_buffer = fopen(outp.fcs_name,"w");  
		fclose(File_buffer);
	}

	if (outp.vel_out)
	{
		strcpy(outp.vel_name,outp.output_files_name);
		strcat (outp.vel_name,".vel");

		File_buffer = fopen(outp.vel_name,"w");
		fclose(File_buffer);
		File_buffer= fopen(outp.vel_name,"a");

		fprintf(File_buffer,"!File contains total atom velocities [h*Angstrom].\n");
		fprintf(File_buffer,"!Total velocities are sum of relative velocity and velocity of center of mass.\n");
		fprintf(File_buffer,"!atom_type atom_number atom_velocity[x] atom_velocity[y] atom_velocity[z]\n");
		fprintf(File_buffer,"!\n");

		fclose(File_buffer);
	}
	if (outp.history_out)
	{
		strcpy(outp.history_name,outp.output_files_name);
		strcat (outp.history_name,".hist");
		File_buffer = fopen(outp.history_name,"w");  
		fclose(File_buffer);

	}
	if (outp.energy_out)
	{
		strcpy(outp.energy_name,outp.output_files_name);
		strcat (outp.energy_name,".eng");
		File_buffer = fopen(outp.energy_name,"w");  
		fclose(File_buffer);
		File_buffer = fopen(outp.energy_name,"a");
		fprintf(File_buffer,"!EKIN EPOT ETOT TKIN\n");
		fclose(File_buffer);
	}
	if (outp.plb_out)
	{
		float plb_header[2];
			  plb_header[0]=(float)(sysp.atom_count_aligned);
		   	  plb_header[1]=-3.0;

		strcpy(outp.plb_name,outp.output_files_name);
		strcat (outp.plb_name,".plb");
		File_buffer = fopen(outp.plb_name,"wb");
		fclose(File_buffer);
		File_buffer= fopen(outp.plb_name,"ab");
		fwrite(plb_header,1,sizeof(plb_header),File_buffer);//write header into .plb file
		fclose(File_buffer);


		strcpy(outp.mol_name,outp.output_files_name);
		strcat (outp.mol_name,".mol");
		File_buffer= fopen(outp.mol_name,"w");
		fclose(File_buffer);
		File_buffer= fopen(outp.mol_name,"a");


		fprintf(File_buffer,"number_of_atoms = %d\n",sysp.atom_count_aligned);
		fprintf(File_buffer,"\n");
		fprintf(File_buffer,"\n");
		fprintf(File_buffer,"atoms\n");
		fprintf(File_buffer,"! i  atom-id  a-type charge chir nbonds bound_atoms\n");

	#if	SUBSTANCE == ARGON

		for (int i=0;i<sysp.atom_count_aligned;i++){
		        fprintf(File_buffer,"%4d %4d-Ar \t Ar      0    0      0\n",i,i);
		}

	#elif SUBSTANCE == SPCE

		for (int i=0;i<sysp.atom_count_aligned;i++){

		        if(i%3==0){
		            fprintf(File_buffer,"%4d %4d-%c \t %c      %9f    0      %d %d %d\n",i,i,'O','O',-0.847600,2,i+1,i+2);
		        }
		        else if (i%3==1){
		            fprintf(File_buffer,"%4d %4d-%c \t %c      %9f    0      %d %d\n",i,i,'H','H',0.423800,1,i-1);
		        }
		        else{
		            fprintf(File_buffer,"%4d %4d-%c \t %c      %9f    0      %d %d\n",i,i,'H','H',0.423800,1,i-2);
		        }
		}

	#elif SUBSTANCE == NITROGEN

		for (int i=0;i<sysp.atom_count_aligned;i++){

			if(i%2==0){
				fprintf(File_buffer,"%4d %4d-%c \t %c      0    0      %d %d\n",i,i,'N','N',1,i+1);
			}
			else{
			    fprintf(File_buffer,"%4d %4d-%c \t %c      0    0      %d %d\n",i,i,'N','N',1,i-1);
			}
		}

	#endif

	fclose(File_buffer);
	}

	if (outp.oob_out)
	{
		strcpy(outp.oob_name,outp.output_files_name);
		strcat (outp.oob_name,".oob");

		File_buffer = fopen(outp.oob_name,"w");
		fclose(File_buffer);
		File_buffer= fopen(outp.oob_name,"a");

		fprintf(File_buffer,"!Box dimensions: lx=%f  ly=%f  lz=%f\n",sysp.lx,sysp.ly,sysp.lz);
		fprintf(File_buffer,"!all atoms contained in outlaying molecule with outlying atoms are printed marked by '*' at the beginning of line\n");
		fprintf(File_buffer,"!Units=angstom\n !density=%f[kg/m^3]\n !number of molecules= %d\n !number of atoms= %d \n",sysp.density,sysp.molecule_count_aligned,sysp.atom_count_aligned);
		fprintf(File_buffer,"!outlying_atom_type outlying_atom_number molecule_number atom[x] atom[y] atom[z] molecule_oficial[x] molecule_oficial[y] molecule_oficial[z] molecule_recalculated[x] molecule_recalculated[y] molecule_recalculated[z]\n");
		fprintf(File_buffer,"!\n");

		fclose(File_buffer);
	}
	if (outp.test4_out)
	{
		strcpy(outp.test4_name,outp.output_files_name);
		strcat (outp.test4_name,".test4");

		File_buffer = fopen(outp.test4_name,"w");
		fclose(File_buffer);
		File_buffer= fopen(outp.test4_name,"a");

		fprintf(File_buffer,"!all atoms contained in outlaying molecule with outlying atoms are printed marked by '*' at the beginning of line\n");
		fprintf(File_buffer,"!Units=angstom\n !density=%f[kg/m^3]\n !number of molecules= %d\n !number of atoms= %d \n",sysp.density,sysp.molecule_count_aligned,sysp.atom_count_aligned);
		fprintf(File_buffer,"!atom_type atom_number molecule_number atom[x] atom[y] atom[z] atom_combiner[x] atom_combiner[y] atom_combiner[z] molecule_oficial[x] molecule_oficial[y] molecule_oficial[z] molecule_recalculated[x] molecule_recalculated[y] molecule_recalculated[z]\n");
		fprintf(File_buffer,"!Mark * atom is out of the box\n");
		fprintf(File_buffer,"!Mark # atom is out of the box with given OVERLAP\n");
		fprintf(File_buffer,"!Mark x relative molecular center of mass is not in zero with given TOLERANCE\n");
		fprintf(File_buffer,"!\n");

		fclose(File_buffer);
	}

	if (outp.fcs_out)
	{
		strcpy(outp.fcs_name,outp.output_files_name);
		strcat(outp.fcs_name,".fcs");
		File_buffer = fopen(outp.fcs_name,"w");  
		fclose(File_buffer);
		File_buffer = fopen(outp.fcs_name,"a");
		fprintf(File_buffer,"!forces acting on each atom, unit [].\n");
		fprintf(File_buffer,"!atom_type atom_number atom_force[x] atom_force[y] atom_force[z]\n");
		fclose(File_buffer);
	}

	if (outp.cp_out)
	{
		strcpy(outp.cp_name,outp.output_files_name);
		strcat (outp.cp_name,".cp");
		File_buffer= fopen(outp.cp_name,"wb");
		fclose(File_buffer);
		File_buffer = fopen(outp.cp_name,"ab");  

		char collumn_3_name[5]="Epot";
		char collumn_4_name[5]="Ekin";
		char collumn_5_name[5]="Pcfg";
		char collumn_6_name[5]="Rho0";
		char collumn_7_name[5]="Ttra";
		char collumn_8_name[5]="Hnc";   //This must be defined for all column except first and second ones
										//max 4 characters here, name of third column, first is Etot and second Tkin by default
		fwrite(&CPmark,1,sizeof(float),File_buffer);
		fwrite(&collumns_in_cp,1,sizeof(int),File_buffer);
		fwrite(&collumn_3_name,sizeof(char),4,File_buffer);
		fwrite(&collumn_4_name,sizeof(char),4,File_buffer);
		fwrite(&collumn_5_name,sizeof(char),4,File_buffer);
		fwrite(&collumn_6_name,sizeof(char),4,File_buffer);
		fwrite(&collumn_7_name,sizeof(char),4,File_buffer);
		fwrite(&collumn_8_name,sizeof(char),4,File_buffer);

	}

	return;
}

void create_log_file(void)
{
	FILE * File_buffer;
	if (outp.log_out)
	{
		strcpy(outp.log_name,outp.output_files_name);
		strcat (outp.log_name,".log");
		File_buffer= fopen(outp.log_name,"w");
		fclose(File_buffer);
		File_buffer= fopen(outp.log_name,"a");
		if (ferror (File_buffer))
		{
			fprintf (stderr,"Error Writing to myfile.txt\n");
			fclose (File_buffer);
		}
		fclose(File_buffer);
	}
	else
	{
		printf(".log file can not be created.");
	}
	return;
}

void set_all_printouts (unsigned int printout_steps_in)
{

	outp.energy_print_period=printout_steps_in;
	outp.history_print_period=printout_steps_in;
	outp.fcs_print_period=printout_steps_in;
	outp.vel_print_period=printout_steps_in;
	outp.oob_print_period=printout_steps_in;
	outp.test4_print_period=printout_steps_in;
	outp.plb_print_period=printout_steps_in;
	outp.cp_print_period=printout_steps_in;
	outp.cfg_print_period=printout_steps_in;
	outp.cfa_print_period=printout_steps_in;

	return;
}
