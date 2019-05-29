/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Xiaojing An <an4@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or process
 * disclosed, or represents that its use would not infringe privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or Lawrence Livermore National Security, LLC. The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <iostream>
#include <assert.h>
#include <utility>
#include <functional>
#include <fstream>      // std::ofstream

#include <ctime>
#include <cstdlib>


// notes for how to setup a good test
// take rank * 100 and make edges between (all local)
// Make one vert per rank a hub.

void usage()  {
    std::cerr << "Usage: -n <int> -f <int> -o <string>\n"
        << " -n <int>    - number of vertices (default 2^17)\n"
        << " -s <int>    - seeds (default time(0))\n"
        << " -f <int>    - range of flags (default 10). Flags generated are: 0, 1, 2, ...\n"
        << " -o <string> - meta data output filename\n"
        << " -h          - print help and exit\n\n";
}

void parse_cmd_line(int argc, char** argv, uint64_t& n, uint64_t& seeds,
        uint64_t& flags, std::string& output_filename) {
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
        std::cout << " " << argv[i];
    }
    std::cout << std::endl;

    bool found_output_filename = false;
    n = (1<<17);
    seeds = time(0);
    flags = 10;

    char c;
    bool prn_help = false;
    while ((c = getopt(argc, argv, "n:s:f:o:h ")) != -1) {
        switch (c) {
            case 'h':  
                prn_help = true;
                break;
            case 'n':
                n = atoll(optarg);
                break;
            case 's':
                seeds = atoll(optarg);
                break;
            case 'f':
                flags = atoll(optarg);
                break; 
            case 'o':
                found_output_filename = true;
                output_filename = optarg;
                break;
            default:
                std::cerr << "Unrecognized option: "<<c<<", ignore."<<std::endl;
                prn_help = true;
                break;
        }
    } 
    if (prn_help || !found_output_filename) {
        usage();
        exit(-1);
    }
}

void uniform_random_flags_generator(uint64_t& n, uint64_t& seeds, uint64_t& flags, std::string& output_filename) {
    std::ofstream ofs (output_filename, std::ofstream::out);

    if (!ofs.is_open()) {
        std::cerr << "Invalid output file" << std::endl;
        exit(-1);
    }

    srand(seeds);
    for(uint64_t i = 0; i < n; i++) {
        ofs << i << " " << (rand() % flags) << std::endl;
    }

    ofs.close();
}

int main(int argc, char** argv) {

    uint64_t      num_vertices;
    uint64_t      seeds;
    uint64_t      flags;
    std::string   output_filename;

    parse_cmd_line(argc, argv, num_vertices, seeds, flags, output_filename);

    std::cout << "Generate Uniformly Random Vertex Matedata"<< std::endl
        << "Output File name = " << output_filename << std::endl;

    uniform_random_flags_generator(num_vertices, seeds, flags, output_filename);

    return 0;
}
