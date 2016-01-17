#include <vector>
#include <libmints/mints.h>

using namespace std;

namespace psi { namespace cctransort {

vector<int> & pitzer2qt(int nirreps, Dimension nmopi, Dimension doccpi, Dimension soccpi, Dimension frzcpi, Dimension frzvpi)
{
  vector<int> order;

  Dimension offset(nirreps);
  offset[0] = 0;
  for(int h=1; h < nirreps; h++)
    offset[h] = offset[h-1] + nmopi[h-1];

  int count = 0; 

  // frozen core orbitals
  for(int h=0; h < nirreps; h++)
    for(int i=0; i < frzcpi[h]; i++)
      order[offset[h]+i] = count++;

  // doubly occupied orbitals
  for(int h=0; h < nirreps; h++)
    for(int i=0; i < (doccpi[h] - frzcpi[h]); i++) // NB: doccpi includes the frozen core orbitals
      order[offset[h]+frzcpi[h]+i] = count++;

  // singly occupied orbitals
  for(int h=0; h < nirreps; h++)
    for(int i=0; i < soccpi[h]; i++)
      order[offset[h]+doccpi[h]+i] = count++;

  // unoccpied orbitals
  for(int h=0; h < nirreps; h++)
    for(int i=0; i < (nmopi[h]-doccpi[h]-soccpi[h]-frzvpi[h]); i++)
      order[offset[h]+doccpi[h]+soccpi[h]+i] = count++;

  // frozen virtual orbitals
  for(int h=0; h < nirreps; h++)
    for(int i=0; i < frzvpi[h]; i++)
      order[offset[h]+nmopi[h]-frzvpi[h]+i] = count++;

  return order;
}

}} 
