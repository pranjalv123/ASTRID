#ifndef ARGS_HPP
#define ARGS_HPP

#include "help.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <cstring>

struct Args {
  std::string infile;
  std::string outfile;
  std::string multindfile;
  std::string cachefile;
  bool cache;

  double constant;
  std::vector<std::string> dms;

  bool octal;

  std::vector<std::string> java_opts;
  bool stacksize_set = false;
  
  Args(int argc, char** argv) {
    constant = 0;
    octal = false;
    cache = false;
    for (int i = 1; i < argc; i++) {
      std::string arg(argv[i]);
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
      else if (arg == "--rapidnj") {
        dms.push_back("rapidnj");
      }
      else if (arg == "--bionj") {
        dms.push_back("bionj");
      }
      else if (arg == "--cache") {
        cachefile = argv[i+1];
        cache = true;
        i++;
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
        std::cerr << help << std::endl;
	exit(1);
      } else if (arg == "-j" || arg == "--java") {
	java_opts.push_back(argv[i+1]);
	if (strncmp(argv[i+1], "-Xss", strlen("-Xss")) == 0) {
	  stacksize_set = true;
	}
	i++;
      }
      else {
        std::cerr << "Unrecognized argument " << arg << std::endl;
	exit(1);
      }


    }
    if (!stacksize_set) {
      java_opts.push_back("-Xss2M");
    }
    if (infile == "") {
      std::cerr << help << std::endl;;
        exit(1);
    }
    if (outfile == "") {
      outfile = infile + ".astrid";
    }

    if (dms.empty()) {
      // by default, we use UPGMA*, then FastME+NNIs, then FastME+SPRs
      dms.push_back("upgma");
      dms.push_back("fastme_nni");
      dms.push_back("fastme_spr");
    }
  }
  
};


#endif
