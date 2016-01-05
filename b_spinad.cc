/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <libdpd/dpd.h>

namespace psi { namespace cctransort {

void b_spinad(void)
{
  dpdbuf4 B, Bs, Ba;

  // This should probably replace the (VV|VV) sort in sort_tei so that we limit disk usage
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
  global_dpd_->buf4_init(&Bs, PSIF_CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
  global_dpd_->buf4_scm(&Bs, 0.0);
  global_dpd_->buf4_init(&Ba, PSIF_CC_BINTS, 0, 9, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
  global_dpd_->buf4_scm(&Ba, 0.0);
  for(int h=0; h < B.params->nirreps; h++) {
    global_dpd_->buf4_mat_irrep_row_init(&B, h);
    global_dpd_->buf4_mat_irrep_row_init(&Bs, h);
    global_dpd_->buf4_mat_irrep_row_init(&Ba, h);
    for(int ab=0; ab < Bs.params->rowtot[h]; ab++) {
      int a = Bs.params->roworb[h][ab][0];
      int b = Bs.params->roworb[h][ab][1];
      int AB = B.params->rowidx[a][b];
      global_dpd_->buf4_mat_irrep_row_rd(&B, h, AB);
      for(int cd=0; cd < Bs.params->coltot[h]; cd++ ) {
        int c = Bs.params->colorb[h][cd][0];
        int d = Bs.params->colorb[h][cd][1];
        int CD = B.params->colidx[c][d];
        int DC = B.params->colidx[d][c];
        Bs.matrix[h][0][cd] = B.matrix[h][0][CD] + B.matrix[h][0][DC];
        Ba.matrix[h][0][cd] = B.matrix[h][0][CD] - B.matrix[h][0][DC];
      }
      global_dpd_->buf4_mat_irrep_row_wrt(&Bs, h, ab);
      global_dpd_->buf4_mat_irrep_row_wrt(&Ba, h, ab);
    }
    global_dpd_->buf4_mat_irrep_row_close(&Ba, h);
    global_dpd_->buf4_mat_irrep_row_close(&Bs, h);
    global_dpd_->buf4_mat_irrep_row_close(&B, h);
  }
  global_dpd_->buf4_close(&Ba);
  global_dpd_->buf4_close(&Bs);
  global_dpd_->buf4_close(&B);
}

}} // namespace psi::cctranssort
