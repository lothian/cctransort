#include <psifiles.h>
#include <libdpd/dpd.h>

namespace psi { namespace cctransort {

double scf_check(int reference)
{
  dpdfile2 f;
  double E1A, E1B, E2AA, E2BB, E2AB;

  if(reference == 2) { // UHF/semicanonical
    global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    E1A = global_dpd_->file2_trace(&f);
    global_dpd_->file2_close(&f);
    global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 2, 2, "fij");
    E1B = global_dpd_->file2_trace(&f);
    global_dpd_->file2_close(&f);

/*
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, "I>J-", "K>L-", "IJ", "KL", 1, "A <IJ|KL>");
    E2AA = global_dpd_->buf4_trace(&A);
    global_dpd_->buf4_close(&A);
*/
    dpdbuf4 ABuf;
    global_dpd_->buf4_init(&ABuf, PSIF_CC_AINTS, 0, "i>j-", "k>l-", "ij", "kl", 1, "A <ij|kl>");
    global_dpd_->buf4_trace(&ABuf);
    global_dpd_->buf4_close(&ABuf);
/*
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, "Ij", "Kl", "Ij", "Kl", 0, "A <Ij|Kl>");
    E2AB = global_dpd_->buf4_trace(&A);
    global_dpd_->buf4_close(&A);
*/

    outfile->Printf("E1A = %20.14f\n", E1A);
    outfile->Printf("E1B = %20.14f\n", E1B);
    outfile->Printf("E2AA = %20.14f\n", E2AA);
    outfile->Printf("E2BB = %20.14f\n", E2BB);
    outfile->Printf("E2AB = %20.14f\n", E2AB);
  }
  else { // RHF/ROHF

  }

  return 0.0;
}

}} // namespace psi::cctransort
