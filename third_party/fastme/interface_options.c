//
//  Copyright 2002-2007 Rick Desper, Olivier Gascuel
//  Copyright 2007-2014 Olivier Gascuel, Stephane Guindon, Vincent Lefort
//
//  This file is part of FastME.
//
//  FastME is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  FastME is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastME.  If not, see <http://www.gnu.org/licenses/>
//


#include "interface_options.h"
#include "config.h"

/*********************************************************/

/*********************************************************/
#ifdef _OPENMP
void Usage (int nbthreads)
#else
void Usage ()
#endif
{
	printf (FLAT"%s\n"BOLD"NAME\n"
		FLAT"\t%s - A comprehensive, accurate and fast distance-based phylogeny inference program.\n\n"
		FLAT"\tVincent Lefort, Richard Desper and Olivier Gascuel,\n"
		FLAT"\tMolecular Biology and Evolution 32(10), 2798-800, 2015.\n"
		FLAT"\tPlease cite this paper if you use this program in your publications.\n", PACKAGE_STRING, PACKAGE_NAME);

	printf (BOLD"\nSYNOPSIS\n"
		BOLD"\t%s"
		FLAT"  ["BOLD"-i "LINE"input data file"FLAT"]"
		FLAT"  ["BOLD"-u "LINE"input user tree file"FLAT"]"
		FLAT"\n\t["BOLD"-o "LINE"output tree file"FLAT"]"
		FLAT"  ["BOLD"-O "LINE"output matrix file"FLAT"]"
		FLAT"  ["BOLD"-I "LINE"output information file"FLAT"]"
		FLAT"\n\t["BOLD"-B "LINE"output bootstrap trees file"FLAT"]"
		FLAT"  ["BOLD"-a"FLAT"]"
		FLAT"\n\t["BOLD"-m "LINE"method"FLAT"]"
		FLAT"  [ "BOLD"-D"LINE"[model]"FLAT" | "BOLD"-P"LINE"[model]"FLAT" ]"
		FLAT"  ["BOLD"-r]"
		FLAT"  ["BOLD"-e]"
		FLAT"  ["BOLD"-g"LINE"[alpha]"FLAT"]"
		FLAT"  ["BOLD"-n"LINE"[NNI]"FLAT"]"
		FLAT"  ["BOLD"-s"FLAT"]"
		FLAT"  ["BOLD"-w "LINE"branch"FLAT"]"
		FLAT"\n\t["BOLD"-d "LINE"datasets"FLAT"]"
		FLAT"  ["BOLD"-b "LINE"replicates"FLAT"]"
		FLAT"  ["BOLD"-z "LINE"seed"FLAT"]"
		FLAT"\n\t["BOLD"-c"FLAT"]"
		FLAT"\n\t["BOLD"-f"FLAT"]"
#ifdef _OPENMP
		FLAT"  ["BOLD"-T "LINE"number of threads"FLAT"]"
#endif
		FLAT"  ["BOLD"-v"FLAT"]"
		FLAT"  ["BOLD"-V"FLAT"]"
		FLAT"  ["BOLD"-h"FLAT"] \n", PACKAGE);

	printf (FLAT"\n\tYou can use %s with no arguments, in this case change the value of\n"
		FLAT"\ta parameter by typing its corresponding character as shown on screen.\n", PACKAGE);

	printf (BOLD"\nOPTIONS\n");

	printf (BOLD"\n\t-i "LINE"input data file"BOLD", --input_data="LINE"input data file"
		FLAT"\n\t\tThe "LINE"input data file"FLAT" contains sequence alignment(s)"
		FLAT"\n\t\tor a distance matrix(ces).\n");

	printf (BOLD"\n\t-u "LINE"input user tree file"BOLD", --user_tree="LINE"input user tree file"
		FLAT"\n\t\t"BOLD"%s "FLAT"may use an existing topology available in the "LINE"input user tree file"
		FLAT"\n\t\twhich corresponds to the input dataset.\n", PACKAGE_NAME);

	printf (BOLD"\n\t-o "LINE"output tree file"BOLD", --output_tree="LINE"output tree file"
		FLAT"\n\t\t"BOLD"%s "FLAT"will write the infered tree into the "LINE"output tree file"FLAT".\n", PACKAGE_NAME);

	printf (BOLD"\n\t-O "LINE"output matrix file"BOLD", --output_matrix="LINE"output matrix file"
		FLAT"\n\t\tUse this option if you want "BOLD"%s "FLAT"to write the distances"
		FLAT"\n\t\tmatrix computed from the input alignment in the "LINE"output matrix file"FLAT".\n", PACKAGE_NAME);

	printf (BOLD"\n\t-I "LINE"output information file"BOLD", --output_info="LINE"output information file"
		FLAT"\n\t\tUse this option if you want "BOLD"%s "FLAT"to write information"
		FLAT"\n\t\tabout its execution in the "LINE"output information file"FLAT".\n", PACKAGE_NAME);

	printf (BOLD"\n\t-B "LINE"output bootstrap trees file"BOLD", --output_boot="LINE"output bootstrap trees file"
		FLAT"\n\t\tUse this option if you want "BOLD"%s "FLAT"to write bootstrap trees"
		FLAT"\n\t\tin the "LINE"bootstrap trees file"FLAT".\n", PACKAGE_NAME);

	printf (BOLD"\n\t-a, --append"
		FLAT"\n\t\tUse this option to append results to existing output files (if any)."
		FLAT"\n\t\tBy default output files will be overwritten.\n");

	printf (BOLD"\n\t-m "LINE"method"BOLD", --method="LINE"method"
		FLAT"\n\t\t"BOLD"%s "FLAT"computes a tree using a distance algorithm."
		FLAT"\n\t\tYou may choose this "LINE"method"FLAT" from:"
		FLAT"\n\t\t"BOLD"TaxAdd_(B)alME"FLAT", "BOLD"TaxAdd_(O)LSME"FLAT", "BOLD"B(I)ONJ"FLAT" (default),"
		FLAT"\n\t\t"BOLD"(N)J"FLAT" or "BOLD"(U)NJ"FLAT".\n", PACKAGE_NAME);

	printf (BOLD"\n\t-d"LINE"[model]"BOLD", --dna="LINE"[model]"
		FLAT"\n\t\tUse this option if your input data file contains DNA sequences alignment(s)."
		FLAT"\n\t\tYou may also indicate the evolutionary "LINE"[model]"FLAT" which can be choosen from:"
		FLAT"\n\t\t"BOLD"(p)-distance"FLAT", "BOLD"R(Y) symmetric"FLAT", "BOLD"(R)Y"FLAT", "BOLD"(J)C69"FLAT", "BOLD"(K)2P"
		FLAT", "BOLD"F8(1)"FLAT", "BOLD"F8(4)"FLAT" (default), "BOLD"(T)N93"FLAT", "BOLD"(L)ogDet"FLAT".\n");

	printf (BOLD"\n\t-p"LINE"[model]"BOLD", --protein="LINE"[model]"
		FLAT"\n\t\tUse this option if your input data file contains protein sequences alignment(s)."
		FLAT"\n\t\tYou may also indicate the evolutionary "LINE"[model]"FLAT" which can be choosen from:"
		FLAT"\n\t\t"BOLD"(p)-distance"FLAT", "BOLD"(F)81 like"FLAT", "BOLD"(L)G"FLAT" (default), "BOLD"(W)AG"FLAT", "BOLD"(J)TT"FLAT", "BOLD"Day(h)off"FLAT", "
		FLAT"\n\t\t"BOLD"(D)CMut"FLAT", "BOLD"(C)pRev"FLAT", "BOLD"(M)tREV"FLAT", "BOLD"(R)tREV"FLAT", "BOLD"HIV(b)"FLAT", "BOLD"H(I)Vw"FLAT" or "BOLD"FL(U)"FLAT".\n");

	printf (BOLD"\n\t-r, --remove_gap"
		FLAT"\n\t\tUse this option to completely remove any site which has a gap in"
		FLAT"\n\t\tany sequence. By default, "BOLD"%s "FLAT"is doing pairwise deletion of gaps.\n", PACKAGE_NAME);

	printf (BOLD"\n\t-e, --equilibrium"
		FLAT"\n\t\tThe equilibrium frequencies for DNA are always estimated by counting"
		FLAT"\n\t\tthe occurence of the nucleotides in the input alignment."
		FLAT"\n\t\tFor amino-acid sequences, the equilibrium frequencies are estimated"
		FLAT"\n\t\tusing the frequencies defined by the substitution model."
		FLAT"\n\t\tUse this option if you whish to estimate the amino-acid frequencies"
		FLAT"\n\t\tby counting their occurence in the input alignment.\n");

	printf (BOLD"\n\t-g"LINE"[alpha]"BOLD", --gamma="LINE"[alpha]"
		FLAT"\n\t\tUse this option if you wish to have gamma distributed rates across sites."
		FLAT"\n\t\tBy default, %s runs with no gamma variation."
		FLAT"\n\t\tIf running %s with gamma distributed rates across sites, the "LINE"[alpha]"FLAT" default value is 1.0."
		FLAT"\n\t\tOnly helpful when the input data file contains sequences alignment(s).\n", PACKAGE_NAME, PACKAGE_NAME);

	printf (BOLD"\n\t-n"LINE"[NNI]"BOLD", --nni="LINE"[NNI]"
		FLAT"\n\t\tUse this option to do "LINE"[NNI]"FLAT" tree topology improvement."
		FLAT"\n\t\tYou may choose the "LINE"[NNI]"FLAT" type from:"
		FLAT"\n\t\t"BOLD"NNI_(B)alME"FLAT" (default) or "BOLD"NNI_(O)LS"FLAT".\n");

	printf (BOLD"\n\t-s, --spr"
		FLAT"\n\t\tUse this option to do "LINE"SPR"FLAT" tree topology improvement.\n");

	printf (BOLD"\n\t-w "LINE"branch"BOLD", --branch_length="LINE"branch"
		FLAT"\n\t\tUse this option to indicate the "LINE"branch"FLAT" length to assign to the tree."
		FLAT"\n\t\tYou may choose the "LINE"branch"FLAT" length from: "
		BOLD"(B)alLS"FLAT" (default), "BOLD"(O)LS"
		FLAT"\n\t\tor "BOLD"(n)one"FLAT". "BOLD"(n)one "FLAT"is only available with BIONJ, NJ or UNJ."
		FLAT"\n\t\tOnly helpful when not improving the tree topology (no NNI nor SPR).\n");

	printf (BOLD"\n\t-D "LINE"datasets"BOLD", --datasets="LINE"datasets"
		FLAT"\n\t\tUse this option to indicate the number of "LINE"datasets"FLAT" in your input"
		FLAT"\n\t\tdata file. Default value is 1.\n");

	printf (BOLD"\n\t-b "LINE"replicates"BOLD", --bootstrap="LINE"replicates"
		FLAT"\n\t\tUse this option to indicate the number of "LINE"replicates"BOLD" %s "FLAT"will"
		FLAT"\n\t\tdo for bootstrapping. Default value is 0."
		FLAT"\n\t\tOnly helpful when the input data file contains sequences alignment(s).\n", PACKAGE_NAME);

	printf (BOLD"\n\t-z "LINE"seed"BOLD", --seed="LINE"seed"
		FLAT"\n\t\tUse this option to initialize randomization with "LINE"seed"FLAT" value."
		FLAT"\n\t\tOnly helpful when bootstrapping.\n");

	printf (BOLD"\n\t-c"
		FLAT"\n\t\tUse this option if you want %s only to compute distance matrix."
		FLAT"\n\t\tOnly helpful when the input data file contains sequences alignment(s).\n", PACKAGE_NAME);

	printf (BOLD"\n\t-f "LINE"number of digits"
		FLAT"\n\t\tUse this option to set the number of digits after the dot to use on output."
		FLAT"\n\t\tDefault precision is 12.\n");

#ifdef _OPENMP
	printf (BOLD"\n\t-T "LINE"number of threads"BOLD", --nb_threads="LINE"number of threads"
		FLAT"\n\t\tUse this option to set the number of threads to use."
		FLAT"\n\t\tDefault "LINE"number of threads"FLAT" is %d.\n", nbthreads);
#endif

	printf (BOLD"\n\t-v "LINE"value"BOLD", --verbose="LINE"value"
		FLAT"\n\t\tSets the verbose level to "LINE"value"FLAT" [0-3]."
		FLAT"\n\t\tDefault "LINE"value"FLAT" is 0.\n");

	printf (BOLD"\n\t-V, --version"
		FLAT"\n\t\tPrints the %s version.\n", PACKAGE_NAME);

	printf (BOLD"\n\t-h, --help"
		FLAT"\n\t\tDisplay this usage.\n\n");

/*
	printf (BOLD"\n\t-t, --TBR"
		FLAT"\n\t\tUse this option to do "LINE"TBR"FLAT" postprocessing.\n");
*/

	return;
}

/*********************************************************/

void Set_Defaults_Input (Options *input)
{
	input->I_data_file	 = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->I_tree_file	 = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->O_tree_file	 = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->O_mat_file	 = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->O_stat_file	 = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->O_boot_file	 = (char *) mCalloc (MAX_FILE_NAME_LENGTH, sizeof(char));
	input->open_mode	 = (char *) mCalloc (3, sizeof(char));

	input->fpI_data_file	= NULL;
	input->fpI_tree_file	= NULL;
	input->fpO_tree_file	= NULL;
	input->fpO_mat_file		= NULL;
	input->fpO_stat_file	= NULL;
	input->fpO_boot_file	= NULL;

	strncpy (input->open_mode, "w", 3);
	input->use_O_mat_file	= FALSE;
	input->is_interleaved	= TRUE;
	input->nb_datasets		= 1;
	input->nb_bootstraps	= 0;
	input->input_type		= MATRIX;
	input->method			= BIONJ;
	input->model			= NONE;
	input->global_aa_fq		= TRUE;
	input->use_gamma		= FALSE;
	input->gamma			= 1.0;
	input->only_mat			= FALSE;
	input->precision		= 12;
	input->use_NNI			= FALSE;
	input->NNI				= BALNNI;
	input->branch			= NONE;
	input->seed				= time (NULL);
	input->no_gap			= FALSE;
	input->use_SPR			= FALSE;
//	input->use_TBR			= FALSE;
	verbose					= 0;

#ifdef _OPENMP
	input->nb_threads		= omp_get_max_threads();
#else
	input->nb_threads		= 1;
#endif

	return;
}

