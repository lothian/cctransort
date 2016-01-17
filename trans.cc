#include <vector>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libtrans/integraltransform.h>
#include <libdpd/dpd.h>

namespace psi { namespace cctransort {

double trans(boost::shared_ptr<Wavefunction> ref, Options& options, std::vector<boost::shared_ptr<MOSpace> > spaces)
{
  IntegralTransform *ints;
  if(options.get_str("REFERENCE") == "RHF")
    ints = new IntegralTransform(ref, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly);
  else if(options.get_str("REFERENCE") == "ROHF")
    ints = new IntegralTransform(ref, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly);
  else if(options.get_str("REFERENCE") == "UHF")
    ints = new IntegralTransform(ref, spaces, IntegralTransform::Unrestricted, IntegralTransform::DPDOnly);
  else
    throw PSIEXCEPTION("Invalid choice of reference wave function.");

  dpd_set_default(ints->get_dpd_id());
  ints->set_keep_dpd_so_ints(true);

  outfile->Printf("\tTransforming integrals...\n");
  outfile->Printf("\t(OO|OO)...\n");
  ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ, IntegralTransform::MakeAndKeep);
  outfile->Printf("\t(OO|OV)...\n");
  ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::vir, IntegralTransform::ReadAndKeep);
  outfile->Printf("\t(OO|VV)...\n");
  ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);

  outfile->Printf("\t(OV|OO)...\n");
  ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::occ, IntegralTransform::MakeAndKeep);
  outfile->Printf("\t(OV|OV)...\n");
  ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir, IntegralTransform::ReadAndKeep);
  outfile->Printf("\t(OV|VV)...\n");
  ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);

  ints->set_keep_dpd_so_ints(false);
  outfile->Printf("\t(VV|OO)...\n");
  ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::occ, IntegralTransform::MakeAndKeep);
  outfile->Printf("\t(VV|OV)...\n");
  ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::vir, IntegralTransform::ReadAndKeep);
  outfile->Printf("\t(VV|VV)...\n");
  ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);

  double efzc = ints->get_frozen_core_energy();

  delete ints;

  return efzc;
}

}} // namespace psi::cctransort
