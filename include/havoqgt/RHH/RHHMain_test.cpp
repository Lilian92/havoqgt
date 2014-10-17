/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
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

#include <iostream>
#include <cstdint>
#include <cassert>
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/allocators/adaptive_pool.hpp>
#include <stdio.h>
#include <time.h>
#include <fstream>

#include "RHHAllocHolder.hpp"
#include "RHHMain.hpp"

using namespace std;

int main (void)
{
  boost::interprocess::file_mapping::remove( "/home/iwabuchi/tmp/test_havoc/tmp/tmp.dat" );
  boost::interprocess::managed_mapped_file mfile( boost::interprocess::create_only, "/home/iwabuchi/tmp/test_havoc/tmp/tmp.dat", 1ULL<<35ULL );
  RHH::AllocatorsHolder holder = RHH::AllocatorsHolder(mfile.get_segment_manager());
  RHH::RHHMain<uint64_t, uint64_t> *rhh = new RHH::RHHMain<uint64_t, uint64_t>(holder, 2ULL);

  srand(1);
  /// srand((unsigned int)time(NULL));

  ofstream input_file;
  input_file.open ("/home/iwabuchi/tmp/test_havoc/tmp/input.txt");

  uint64_t num_elems = 0;
  for (uint64_t i = 0; i < (1<<24ULL); ++i) {
    if (i % (1<<20ULL) == 0) std::cout << "!";
    uint64_t key = rand() % (1<<20ULL);
    uint64_t val = rand() % (1<<16ULL);
    bool result = rhh->insert_uniquely(holder, key, val);
    //std::cout << key << "\t" << val << "\t" << result << std::endl;
    input_file << key << "\t" << val << "\t" << result << std::endl;
    num_elems += result;
  }
  input_file.close();

  std::cout << "\n=====================================\n";
  ofstream output_file;
  output_file.open ("/home/iwabuchi/tmp/test_havoc/tmp/output.txt");
  rhh->disp_elems(output_file);
  output_file.close();

  std::cout << "Segment size = " << mfile.get_segment_manager()->get_size() << std::endl;
  std::cout << "Used size = " << mfile.get_segment_manager()->get_size() - mfile.get_segment_manager()->get_free_memory() << std::endl;

  delete rhh;

}
