#ifndef ARGS_HPP
#define ARGS_HPP

#include "help.hpp"

struct Args {
  string infile;
  string outfile;
  string multindfile;
  
  double constant;
  vector<string> dms;

  bool octal;

  Args(int argc, char** argv) {
    constant = 0;
    octal = false;
    for (int i = 1; i < argc; i++) {
      string arg(argv[i]);
      if (arg == "-i" || arg == "--input") {
        infile = argv[i+1];
        i++;
      }

      else if (arg == "-o" || arg == "--output") {
        outfile = argv[i+1];
        i++;
      }

      else if (arg == "-u") {
        dms.push_back("upgma");
      }

      else if (arg == "-f") {
        dms.push_back("fastme");
      }

      else if (arg == "-n") {
        dms.push_back("fastme_nni");
      }

      else if (arg == "-s") {
        dms.push_back("fastme_spr");
      }
      else if (arg == "--bionj") {
        dms.push_back("bionj");
      }

      else if (arg == "-x") {
        constant = atof(argv[i+1]);
        i++;
      }

      else if (arg == "-c") {
        octal=true;
      }
      else if (arg == "-a" || arg == "--multind") {
        multindfile = argv[i+1];
	i++;
      }
      else if (arg == "-h" || arg == "--help") {
	cerr << help << endl;
	exit(1);
      }
      else {
	cerr << "Unrecognized argument " << arg << endl;
	exit(1);
      }


    }
    if (infile == "") {
      cerr << help << endl;;
        exit(1);
    }
    if (outfile == "") {
      outfile = infile + ".astrid";
    }

    if (dms.empty()) {
      dms.push_back("auto");
    }
  }
  
};


#endif
