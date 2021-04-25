#include<cstdlib>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include <algorithm>
#include <math.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<iomanip>
#include<unistd.h>

using namespace std;


//sequence
string fasta;
//sequence route
string fasta_route;
//seqence length
int fasta_length;
//contact map route
string dmap_route;

//fragment library route
string frag_lib_3_route;
string frag_lib_9_route; 

//rosetta executable file
string rosetta_route;

//check whether fastaFile is input
bool fFile;
//check whether cmapFile is input
bool dFile;
//check whether frag_library are input
bool f3File;
bool f9File;
//check whether rosetta executable file is input
bool rFile;


void getParameter(int argc, char ** argv, int & i);
void extractParameters(int argc, char ** argv);
bool is_amino(char &i);
// read the fasta file
void getFasta();

void apply();


int main(int argc, char ** argv){
	cout << endl;
	cout << "###########################################################################" << endl;
	cout << "#                                CGLFold                                  #" << endl;
	cout << "#                              Version: 1.0                               #" << endl;
	cout << "#          A contact-assisted de novo protein structure prediction        #" << endl;
	cout << "#    using global exploration and loop perturbation sampling algorithm    #" << endl;
	cout << "#                         Copyright (C) Guijun Zhang                      #" << endl;
	cout << "#                     College of Information engineering                  #" << endl;
	cout << "#            Zhejiang university of technology, Hangzhou, China           #" << endl;
	cout << "###########################################################################" << endl;
	

	extractParameters(argc, argv);
	
	getFasta();
	
	apply();
	
	cout << endl;
	cout << "###########################################################################" << endl;
	cout << "#                                  END                                    #" << endl;
	cout << "###########################################################################" << endl;
  
	return 0;
}

void getFasta(){
    ifstream fasta_file( fasta_route.c_str() );
    
    if ( fasta_file.is_open() ){
	string line;
	while( getline(fasta_file, line) ){
	    if ( line[0] != '>' )
		fasta += line;
	}
	fasta_file.close();
	
	for (int i = 0; i < fasta.size(); ++i){
	    if ( fasta[i] == ' ' ){
		    fasta.erase(i,1);
		    --i;
	    }
	    else{
		if ( !is_amino( fasta[i] ) ){
		    cout << "Error! Invalid amino acid " << fasta[i] << endl;
		    exit(0);
		}
	    }
	}
	fasta_length = fasta.size();
	cout << "#  " << fasta << endl;
	cout << "#  Fasta length: " << fasta_length << endl;
    }
    else{
	cout << "Error! fasta file can not open " << fasta_route << endl;
	exit(0);
    }
}

bool is_amino(char &r){
	if (r == 'A')
            return true;
        else if (r == 'C')
            return true;
        else if (r == 'D')
            return true;
        else if (r == 'E')
            return true;
        else if (r == 'F')
            return true;
        else if (r == 'G')
            return true;
        else if (r == 'H')
            return true;
        else if (r== 'I')
            return true;
        else if (r == 'K')
            return true;
        else if (r == 'L')
            return true;
        else if (r == 'M')
            return true;
        else if (r == 'N')
            return true;
        else if (r == 'P')
            return true;
        else if (r == 'Q')
            return true;
        else if (r == 'R')
            return true;
        else if (r == 'S')
            return true;
        else if (r == 'T')
            return true;
        else if (r == 'V')
            return true;
        else if (r == 'W')
            return true;
        else if (r == 'Y')
            return true;
        else 
	    return false;
}

void extractParameters(int argc, char ** argv) {
    int i = 1;
    while (i < argc)
        getParameter(argc, argv, i);
    if (!fFile) {
        cout << endl;
        cout << "#  Error! fasta File must be provided" << endl << endl;
        exit(0);
    }
    if (!dFile) {
        cout << endl;
        cout << "#  Error! inter-residue distance File must be provided" << endl << endl;
        exit(0);
    }
    if (!f3File) {
        cout << endl;
        cout << "#  Error! 3-mer fragment library File must be provided" << endl << endl;
        exit(0);
    }
    if (!f9File) {
        cout << endl;
        cout << "#  Error! 9-mer fragment library File must be provided" << endl << endl;
        exit(0);
    }
    if (!rFile) {
        cout << endl;
        cout << "#  Error! rosetta executable File must be provided" << endl << endl;
        exit(0);
    }
}

void getParameter(int argc, char ** argv, int & i){
	string flag( argv[i] );
	if ( flag == "-f" ){
		if (argc < i + 2){
			cout << endl;
			cout << "#  Error! No fasta file provided" << endl << endl;
			cout << "   -f 		fasta  			: fasta file  " << endl;
			cout << "   -d 		distance   		: inter-residue distance file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			cout << "   -rosetta	rosetta executable file : executable file of ClassicAbinitio protocol of Rosetta" << endl;
			exit(0);
		}
		fasta_route = argv[++i];
		fFile = true;
	}
	else if( flag == "-d" ){
		if (argc < i + 2){
			cout << endl;
			cout << "#  Error! No fasta file provided" << endl << endl;
			cout << "   -f 		fasta  			: fasta file  " << endl;
			cout << "   -d 		distance   		: inter-residue distance file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			cout << "   -rosetta	rosetta executable file : executable file of ClassicAbinitio protocol of Rosetta" << endl;
			exit(0);
		}
		dmap_route = argv[++i];
		dFile = true;
	}
	else if( flag == "-frag3" ){
		if (argc < i + 2){
			cout << endl;
			cout << "#  Error! No fasta file provided" << endl << endl;
			cout << "   -f 		fasta  			: fasta file  " << endl;
			cout << "   -d 		distance   		: inter-residue distance file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			cout << "   -rosetta	rosetta executable file : executable file of ClassicAbinitio protocol of Rosetta" << endl;
			exit(0);
		}
		frag_lib_3_route = argv[++i];
		f3File = true;
	}
	else if( flag == "-frag9" ){
		if (argc < i + 2){
			cout << endl;
			cout << "#  Error! No fasta file provided" << endl << endl;
			cout << "   -f 		fasta  			: fasta file  " << endl;
			cout << "   -d 		distance   		: inter-residue distance file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			cout << "   -rosetta	rosetta executable file : executable file of ClassicAbinitio protocol of Rosetta" << endl;
			exit(0);
		}
		frag_lib_9_route = argv[++i];
		f9File = true;
	}
	else if( flag == "-rosetta" ){
		if (argc < i + 2){
			cout << endl;
			cout << "#  Error! No rosetta software provided" << endl << endl;
			cout << "   -f 		fasta  			: fasta file  " << endl;
			cout << "   -d 		distance   		: inter-residue distance file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			cout << "   -rosetta	rosetta executable file : executable file of ClassicAbinitio protocol of Rosetta" << endl;
			exit(0);
		}
		rosetta_route = argv[++i];
		rFile = true;
	}
	++i;
}

void apply(){
	stringstream ddd;
	ddd << "cp " << dmap_route.c_str() << " ./dist.txt";
	string ddd_2;
	getline(ddd, ddd_2);
	system( ddd_2.c_str() );
	stringstream command;
	command << rosetta_route << " "
		<< "-in:file:fasta " << fasta_route << " "
		<< "-in:file:frag3 " << frag_lib_3_route << " "
		<< "-in:file:frag9 " << frag_lib_9_route << " "
		<< "-abinitio:increase_cycles " << 1 << " "
		<< "-nstruct " << 1 << " "
		<< "-out:pdb";
	string run_rosetta;
	getline(command, run_rosetta);
		
	system(run_rosetta.c_str());
	
	system("rm S_00000001.pdb score.fsc default.out dist.txt");
}

