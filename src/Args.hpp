#ifndef ARGS_HPP
#define ARGS_HPP

struct Args {
  string infile;
  string outfile;

  double constant;
  vector<string> dms;

  bool octal, expand;

  Args(int argc, char** argv) {
    constant = 0;
    octal = false;
    expand = false;
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

      else if (arg == "-x") {
        constant = atof(argv[i+1]);
        i++;
      }

      else if (arg == "-c") {
        octal=true;
      }

      else if (arg == "-e") {
        expand=true;
      }


    }
    if (infile == "") {
        cerr << "Input file required";
        exit(1);
    }
    if (outfile == "") {
      outfile = infile + ".astrid";
    }

  }
};


#endif
