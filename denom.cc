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

void denom_rhf(vector(int> &openpi);
void denom_uhf(vector(int> &openpi);

void denom(int reference, vector<int> &openpi)
{
  if(reference == 2) denom_uhf(vector<int> &openpi);
  else denom_rhf(vector<int> &openpi);
}

void denom_uhf(vector<int> &openpi)
{
  dpdfile2 d1, fIJ, fij, fAB, fab;
  dpdbuf4 d2;

  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_mat_init(&fIJ);
  global_dpd_->file2_mat_rd(&fIJ);
  
  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
  global_dpd_->file2_mat_init(&fij);
  global_dpd_->file2_mat_rd(&fij);
  
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_mat_init(&fAB);
  global_dpd_->file2_mat_rd(&fAB);
  
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
  global_dpd_->file2_mat_init(&fab);
  global_dpd_->file2_mat_rd(&fab);

  global_dpd_->file2_init(&d1, PSIF_CC_OEI, 0, 0, 1, "dIA");
  global_dpd_->file2_mat_init(&d1);
  for(int h=0; h < d1.params->nirreps; h++) {
      for(int i=0; i < d1.params->ppi[h]; i++) {
	  double fii = fIJ.matrix[h][i][i];
	  for(int a=0; a < d1.params->qpi[h]; a++) {
	      double faa = fAB.matrix[h][a][a];
	      d1.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
    }
  global_dpd_->file2_mat_wrt(&d1);
  global_dpd_->file2_mat_close(&d1);
  global_dpd_->file2_close(&d1);

  global_dpd_->file2_init(&d1, PSIF_CC_OEI, 0, 2, 3, "dia");
  global_dpd_->file2_mat_init(&d1);
  for(int h=0; h < d1.params->nirreps; h++) {
      for(int i=0; i < d1.params->ppi[h]; i++) {
	  double fii = fij.matrix[h][i][i];
	  for(int a=0; a < d1.params->qpi[h]; a++) {
	      double faa = fab.matrix[h][a][a];
	      d1.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
    }
  global_dpd_->file2_mat_wrt(&d1);
  global_dpd_->file2_mat_close(&d1);
  global_dpd_->file2_close(&d1);

  global_dpd_->buf4_init(&d2, PSIF_CC_DENOM, 0, "I>J+", "A>B+", 0, "dIJAB");
  for(int h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&d2, h);
      for(int ij=0; ij < d2.params->rowtot[h]; ij++) {
	  int i = d2.params->roworb[h][ij][0];
	  int j = d2.params->roworb[h][ij][1];
	  int isym = d2.params->psym[i];
	  int jsym = d2.params->qsym[j];
	  int I = i - d2.params->poff[isym];
	  int J = j - d2.params->qoff[jsym];
	  double fii = fIJ.matrix[isym][I][I];
	  double fjj = fIJ.matrix[jsym][J][J];

	  for(int ab=0; ab < d2.params->coltot[h]; ab++) {
	      int a = d2.params->colorb[h][ab][0];
	      int b = d2.params->colorb[h][ab][1];
	      int asym = d2.params->rsym[a];
	      int bsym = d2.params->ssym[b];
	      int A = a - d2.params->roff[asym];
	      int B = b - d2.params->roff[bsym];
	      double faa = fAB.matrix[asym][A][A];
	      double fbb = fAB.matrix[bsym][B][B];

	      d2.matrix[h][ij][ab] = 1.0/(fii + fjj - faa - fbb);
	    }
	}
      global_dpd_->buf4_mat_irrep_wrt(&d2, h);
      global_dpd_->buf4_mat_irrep_close(&d2, h);
    }
  global_dpd_->buf4_close(&d2);

  global_dpd_->buf4_init(&d2, PSIF_CC_DENOM, 0, "i>j+", "a>b+", 0, "dijab");
  for(int h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&d2, h);
      for(int ij=0; ij < d2.params->rowtot[h]; ij++) {
          int i = d2.params->roworb[h][ij][0];
          int j = d2.params->roworb[h][ij][1];
          int isym = d2.params->psym[i];
          int jsym = d2.params->qsym[j];
          int I = i - d2.params->poff[isym];
          int J = j - d2.params->qoff[jsym];
          double fii = fij.matrix[isym][I][I];
          double fjj = fij.matrix[jsym][J][J];

          for(int ab=0; ab < d2.params->coltot[h]; ab++) {
              int a = d2.params->colorb[h][ab][0];
              int b = d2.params->colorb[h][ab][1];
              int asym = d2.params->rsym[a];
              int bsym = d2.params->ssym[b];
              int A = a - d2.params->roff[asym];
              int B = b - d2.params->roff[bsym];
              double faa = fab.matrix[asym][A][A];
              double fbb = fab.matrix[bsym][B][B];

              d2.matrix[h][ij][ab] = 1.0/(fii + fjj - faa - fbb);
            }
        }
      global_dpd_->buf4_mat_irrep_wrt(&d2, h);
      global_dpd_->buf4_mat_irrep_close(&d2, h);
    }
  global_dpd_->buf4_close(&d2);

  global_dpd_->buf4_init(&d2, PSIF_CC_DENOM, 0, "Ij", "Ab", 0, "dIjAb");
  for(int h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&d2, h);
      for(int ij=0; ij < d2.params->rowtot[h]; ij++) {
          int i = d2.params->roworb[h][ij][0];
          int j = d2.params->roworb[h][ij][1];
          int isym = d2.params->psym[i];
          int jsym = d2.params->qsym[j];
          int I = i - d2.params->poff[isym];
          int J = j - d2.params->qoff[jsym];
          double fii = fIJ.matrix[isym][I][I];
          double fjj = fij.matrix[jsym][J][J];

          for(int ab=0; ab < d2.params->coltot[h]; ab++) {
              int a = d2.params->colorb[h][ab][0];
              int b = d2.params->colorb[h][ab][1];
              int asym = d2.params->rsym[a];
              int bsym = d2.params->ssym[b];
              int A = a - d2.params->roff[asym];
              int B = b - d2.params->roff[bsym];
              double faa = fAB.matrix[asym][A][A];
              double fbb = fab.matrix[bsym][B][B];

              d2.matrix[h][ij][ab] = 1.0/(fii + fjj - faa - fbb);
            }
        }
      global_dpd_->buf4_mat_irrep_wrt(&d2, h);
      global_dpd_->buf4_mat_irrep_close(&d2, h);
    }
  global_dpd_->buf4_close(&d2);

  global_dpd_->file2_mat_close(&fIJ);
  global_dpd_->file2_mat_close(&fij);
  global_dpd_->file2_mat_close(&fAB);
  global_dpd_->file2_mat_close(&fab);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fij);
  global_dpd_->file2_close(&fAB);
  global_dpd_->file2_close(&fab);
}

void denom_rhf(vector<int> &openpi)
{
  dpdfile2 fIJ, fij, fAB, fab;
  dpdfile2 d1;
  dpdbuf4 d2;

  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_mat_init(&fIJ);
  global_dpd_->file2_mat_rd(&fIJ);
  
  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
  global_dpd_->file2_mat_init(&fij);
  global_dpd_->file2_mat_rd(&fij);
  
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_mat_init(&fAB);
  global_dpd_->file2_mat_rd(&fAB);
  
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");
  global_dpd_->file2_mat_init(&fab);
  global_dpd_->file2_mat_rd(&fab);

  global_dpd_->file2_init(&d1, PSIF_CC_OEI, 0, 0, 1, "dIA");
  global_dpd_->file2_mat_init(&dIA);
  for(int h=0; h < d1.params->nirreps; h++) {
      for(int i=0; i < d1.params->ppi[h]; i++) {
	  double fii = fIJ.matrix[h][i][i];
	  for(int a=0; a < d1.params->qpi[h] - openpi[h]; a++) {
	      double faa = fAB.matrix[h][a][a];
	      d1.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
    }
  global_dpd_->file2_mat_wrt(&d1);
  global_dpd_->file2_mat_close(&d1);
  global_dpd_->file2_close(&d1);

  global_dpd_->file2_init(&d1, PSIF_CC_OEI, 0, 0, 1, "dia");
  global_dpd_->file2_mat_init(&d1);
  for(int h=0; h < nirreps; h++) {
      for(int i=0; i < d1.params->ppi[h] - openpi[h]; i++) {
	  double fii = fij.matrix[h][i][i];
	  for(int a=0; a < d1.params->qpi[h]; a++) {
	      double faa = fab.matrix[h][a][a];
	      d1.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
    }
  global_dpd_->file2_mat_wrt(&d1);
  global_dpd_->file2_mat_close(&d1);
  global_dpd_->file2_close(&d1);

  global_dpd_->buf4_init(&d2, PSIF_CC_DENOM, 0, "i>j+", "a>b+", "dIJAB");
  for(int h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&d2, h);
      for(int ij=0; ij < d2.params->rowtot[h]; ij++) {
	  int i = d2.params->roworb[h][ij][0];
	  int j = d2.params->roworb[h][ij][1];
	  int isym = d2.params->psym[i];
	  int jsym = d2.params->qsym[j];
	  int I = i - d2.params->poff[isym];
	  int J = j - d2.params->qoff[jsym];
	  double fii = fIJ.matrix[isym][I][I];
	  double fjj = fIJ.matrix[jsym][J][J];
	  for(int ab=0; ab < d2.params->coltot[h]; ab++) {
	      int a = d2.params->colorb[h][ab][0];
	      int b = d2.params->colorb[h][ab][1];
	      int asym = d2.params->rsym[a];
	      int bsym = d2.params->ssym[b];
	      int A = a - d2.params->roff[asym];
	      int B = b - d2.params->roff[bsym];
	      double faa = fAB.matrix[asym][A][A];
	      double fbb = fAB.matrix[bsym][B][B];
	      d2.matrix[h][ij][ab] = ((A >= (virtpi[asym] - openpi[asym])) || (B >= (virtpi[bsym] - openpi[bsym])) ?
		 0.0 : 1.0/(fii + fjj - faa - fbb));
	    }
	}
      global_dpd_->buf4_mat_irrep_wrt(&d2, h);
      global_dpd_->buf4_mat_irrep_close(&d2, h);
    }
  global_dpd_->buf4_close(&d2);

  global_dpd_->buf4_init(&d2, PSIF_CC_DENOM, 0, "i>j+", "a>b+", 0, "dijab");
  for(int h=0; h < d2.params->nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&d2, h);
      for(int ij=0; ij < d2.params->rowtot[h]; ij++) {
	  int i = d2.params->roworb[h][ij][0];
	  int j = d2.params->roworb[h][ij][1];
	  int isym = d2.params->psym[i];
	  int jsym = d2.params->qsym[j];
	  int I = i - d2.params->poff[isym];
	  int J = j - d2.params->qoff[jsym];
	  double fii = fij.matrix[isym][I][I];
	  double fjj = fij.matrix[jsym][J][J];
	  for(int ab=0; ab < d2.params->coltot[h]; ab++) {
	      int a = d2.params->colorb[h][ab][0];
	      int b = d2.params->colorb[h][ab][1];
	      int asym = d2.params->rsym[a];
	      int bsym = d2.params->ssym[b];
	      int A = a - d2.params->roff[asym];
	      int B = b - d2.params->soff[bsym];
	      double faa = fab.matrix[asym][A][A];
	      double fbb = fab.matrix[bsym][B][B];
	      d2.matrix[h][ij][ab] = ((I >= (occpi[isym] - openpi[isym])) || (J >= (occpi[jsym] - openpi[jsym])) ?
		 0.0 : 1.0/(fii + fjj - faa - fbb));
	    }
	}
      global_dpd_->buf4_mat_irrep_wrt(&d2, h);
      global_dpd_->buf4_mat_irrep_close(&d2, h);
    }
  global_dpd_->buf4_close(&d2);

  global_dpd_->buf4_init(&d2, PSIF_CC_DENOM, 0, "ij", "ab", 0, "dIjAb");
  for(int h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&d2, h);
      for(int ij=0; ij < d2.params->rowtot[h]; ij++) {
	  int i = d2.params->roworb[h][ij][0];
	  int j = d2.params->roworb[h][ij][1];
	  int isym = d2.params->psym[i];
	  int jsym = d2.params->qsym[j];
	  int I = i - d2.params->poff[isym];
	  int J = j - d2.params->qoff[jsym];
	  double fii = fIJ.matrix[isym][I][I];
	  double fjj = fij.matrix[jsym][J][J];
	  for(int ab=0; ab < d2.params->coltot[h]; ab++) {
	      int a = d2.params->colorb[h][ab][0];
	      int b = d2.params->colorb[h][ab][1];
	      int asym = d2.params->rsym[a];
	      int bsym = d2.params->ssym[b];
	      int A = a - d2.params->roff[asym];
	      int B = b - d2.params->soff[bsym];
	      double faa = fAB.matrix[asym][A][A];
	      double fbb = fab.matrix[bsym][B][B];
	      d2.matrix[h][ij][ab] = ((A >= (virtpi[asym] - openpi[asym])) || (J >= (occpi[jsym] - openpi[jsym])) ?
		 0.0 : 1.0/(fii + fjj - faa - fbb));
	    }
	}
      global_dpd_->buf4_mat_irrep_wrt(&d2, h);
      global_dpd_->buf4_mat_irrep_close(&d2, h);
    }
  global_dpd_->buf4_close(&d2);

  global_dpd_->file2_mat_close(&fIJ);
  global_dpd_->file2_mat_close(&fij);
  global_dpd_->file2_mat_close(&fAB);
  global_dpd_->file2_mat_close(&fab);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fij);
  global_dpd_->file2_close(&fAB);
  global_dpd_->file2_close(&fab);

}

}} // namespace psi::ccsort