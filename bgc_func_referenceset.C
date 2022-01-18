#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>
#include <float.h>

#include <math.h>

#include "bgc.h"

using namespace std;

/*------------------------------------------------------*/
/*         data reading functions for reference loci    */
/*------------------------------------------------------*/

/* read in genotype likelihood data for the admixed and parental populations */
void getDataref(string hybridfile, string p0file, string p1file, datacont * data,
		datadimcont * datadim){

  // read and store hybrid data
  readNgsFile(hybridfile,data->refglh,datadim->nrefloci,datadim->nind);
  
  // read and store p0 data
  readNgsFile(p0file,data->refglp0,datadim->nrefloci,datadim->nindP0);
  
  // read and store p1 data
  readNgsFile(p1file,data->refglp1,datadim->nrefloci,datadim->nindP1);
}  

/* get key dimensions, nloci, ninds for gl format ngs data */
void getRefDims(string hybridfile, datadimcont * datadim){

  string line, element;
  ifstream infile;
  istringstream stream;
  
  // hybrids
  infile.open(hybridfile.c_str());
  if (!infile){
    cerr << "Cannot open file " << hybridfile << endl;
    exit(1);
  }
  // read line with data dimensions
  getline(infile, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of  loci
  datadim->nrefloci = atoi(element.c_str());  
  infile.close();

}


/* gibbs sample update of parental population counts */
void updatePcntRef(datacont * data, auxcont * auxvar, datadimcont * datadim, gsl_matrix_int * pcnt, 
		mat_container * gl, int nind, gsl_matrix * pi){
  int i, k, n, a, b;
  double prob[GEN];
  unsigned int geno[GEN];

  // set pcnt to zero
  gsl_matrix_int_set_zero(pcnt);

  for(i=0; i<datadim->nrefloci; i++){
    //cerr << "we have " << nind << " for snps " << i << endl;
    for(n=0; n<nind; n++){
      for(k=0; k<GEN; k++){
	prob[k] = gsl_matrix_get(gl[i].mat, n, k); // store genotype likelihoods
      }
      // modify by allele frequency priors
      prob[0] = prob[0] * gsl_pow_2(gsl_matrix_get(pi, i, 0));
      prob[2] = prob[2] * gsl_pow_2(gsl_matrix_get(pi, i, 1));
      prob[1] = prob[1] * 2 * gsl_matrix_get(pi, i, 0) * gsl_matrix_get(pi, i, 1);

      // sample genotype
      gsl_ran_multinomial(r, GEN, 1, prob, geno);
      for(k=0; k<GEN; k++){
	if(geno[k] == 1){
	  a = gsl_matrix_int_get(pcnt, i, 0);
	  b = gsl_matrix_int_get(pcnt, i, 1);	  
	  gsl_matrix_int_set(pcnt, i, 0, a+k);
	  gsl_matrix_int_set(pcnt, i, 1, b+2-k); // two columns, one for each allele
	}
      }
    }
   // cerr << gsl_matrix_int_get(pcnt, i, 0) << ", ";
  }
  //cerr << endl;
}

/* function to update z, ancestry of hybrids, uses gibbs sampling similar to structure, but with ngs data */
/* we analytically sum across the unkown genotype */
void updateZngsRef(paramcont * param, datacont * data, datadimcont * datadim, auxcont * auxvar){
  int i, n, k, l;
  double p00 = 0, p01 = 0, p10 = 0, p11 = 0, ptot = 0;
  double pr[4]; // array of probabilities
  unsigned int sam[4]; // sampled ancestry for a locus
  int finite;
  double probtemp;

  // loop through loci
  for(i=0; i<datadim->nrefloci; i++){

    // loop through individuals
    for(n=0; n<datadim->nind; n++){
      p00 = 0;
      p01 = 0;
      p11 = 0;
      // first homozygous genotype
      k = 0; l = 0;
      // probability of data given the genotype, i.e., p(d | g) or gl
      probtemp = gsl_matrix_get(data->refglh[i].mat, n, 0);
      // multiply by probability of genotype given ancestry
      // p(d | g) * p(g | z), note k = l
      p00 += probtemp * gsl_matrix_get(param->refpi0Cur, i, k) * gsl_matrix_get(param->refpi0Cur, i, l);
      p01 += probtemp * gsl_matrix_get(param->refpi0Cur, i, k) * gsl_matrix_get(param->refpi1Cur, i, l);
      p11 += probtemp * gsl_matrix_get(param->refpi1Cur, i, k) * gsl_matrix_get(param->refpi1Cur, i, l);

      //cerr << i << "," << n << ": " << p00 << "," << p01 << "," << p11 << endl;
      
      // second homozygous genotype
      k = 1; l = 1;
      // probability of data given the genotype, i.e., p(d | g) or gl
      probtemp = gsl_matrix_get(data->refglh[i].mat, n, 2);
      // multiply by probability of genotype given ancestry
      // p(d | g) * p(g | z), note k = l
      p00 += probtemp * gsl_matrix_get(param->refpi0Cur, i, k) * gsl_matrix_get(param->refpi0Cur, i, l);
      p01 += probtemp * gsl_matrix_get(param->refpi0Cur, i, k) * gsl_matrix_get(param->refpi1Cur, i, l);
      p11 += probtemp * gsl_matrix_get(param->refpi1Cur, i, k) * gsl_matrix_get(param->refpi1Cur, i, l);

      //cerr << i << "," << n << ": " << p00 << "," << p01 << "," << p11 << endl;

      
      // heterozgyous genotype
      k = 0; l = 1;
      // probability of data given the genotype, i.e., p(d | g) or gl
      probtemp = gsl_matrix_get(data->refglh[i].mat, n, 1);
      p00 += probtemp * 2 * gsl_matrix_get(param->refpi0Cur, i, k) * gsl_matrix_get(param->refpi0Cur, i, l);
      p01 += probtemp * (gsl_matrix_get(param->refpi0Cur, i, k) * gsl_matrix_get(param->refpi1Cur, i, l) +
      			 gsl_matrix_get(param->refpi1Cur, i, k) * gsl_matrix_get(param->refpi0Cur, i, l));
      p11 += probtemp * 2 * gsl_matrix_get(param->refpi1Cur, i, k) * gsl_matrix_get(param->refpi1Cur, i, l);

      //cerr << i << "," << n << ": " << p00 << "," << p01 << "," << p11 << endl;
      
      // calculate p(z|phi) and mulitply by p(x|g) p(g|z,pi)
      p00 *= gsl_pow_2(1 - gsl_matrix_get(param->refphi, i, n));
      p01 *= (1 - gsl_matrix_get(param->refphi, i, n)) * gsl_matrix_get(param->refphi, i, n);
      p11 *= gsl_pow_2(gsl_matrix_get(param->refphi, i, n));

      //cerr << i << "," << n << ": " << p00 << "," << p01 << "," << p11 << endl;

      
      //cerr << "ancestry: " << i << " " << n << " " << p00 << "," << p01 << "," << p11 << endl;
      finite = gsl_isinf(p00);
      if (finite == 1){
	p00 = DBL_MAX;
      }
      else if (finite == -1){
	p00 = DBL_MIN;
      }      
      finite = gsl_isinf(p01);
      if (finite == 1){
	p01 = DBL_MAX;
      }
      else if (finite == -1){
	p01 = DBL_MIN;
      }      
      finite = gsl_isinf(p11);
      if (finite == 1){
	p11 = DBL_MAX;
      }
      else if (finite == -1){
	p11 = DBL_MIN;
      }

      p10 = p01;
      ptot = p00 + p01 + p10 + p11;
      //cerr << ptot << endl;
      pr[0] = p00 / ptot;
      pr[1] = p01 / ptot;
      pr[2] = p10 / ptot;
      pr[3] = p11 / ptot;
      //cerr << i << "," << n << ": " << pr[0] << "," << pr[1] << "," << pr[2] << "," << pr[3] << endl << endl;
      // sample ancestry pair based on probabilities
      gsl_ran_multinomial(r, 4, 1, pr, sam);
      //cerr << sam[0] << " " << sam[1] << " " << sam[2] << " " << sam[3] << endl;
      if(sam[0] == 1){ // homo pop 0
	gsl_matrix_int_set(param->refzCur[i].mat, n, 0, 0);
	gsl_matrix_int_set(param->refzCur[i].mat, n, 1, 0);
      }
      else if(sam[1] == 1){ // het
	gsl_matrix_int_set(param->refzCur[i].mat, n, 0, 0);
	gsl_matrix_int_set(param->refzCur[i].mat, n, 1, 1);
      }
      else if(sam[2] == 1){ // het
	gsl_matrix_int_set(param->refzCur[i].mat, n, 0, 1);
	gsl_matrix_int_set(param->refzCur[i].mat, n, 1, 0);
      }
      else{ // homo pop 1
	gsl_matrix_int_set(param->refzCur[i].mat, n, 0, 1);
	gsl_matrix_int_set(param->refzCur[i].mat, n, 1, 1);
      }
    }
  }
}

/* function to update hi, uses gibbs sampling and stores results in a hist */
void updateHiRef(paramcont * param, datadimcont * datadim, int allDiploid){
  int n, i;
  double a, b, h;
  double prior = 0.5;
  
  // loop through individuals
  for(n=0; n<datadim->nind; n++){
    a = 0;
    for(i=0; i<datadim->nrefloci; i++){
      a += gsl_matrix_int_get(param->refzCur[i].mat, n, 0) +
	gsl_matrix_int_get(param->refzCur[i].mat, n, 1);
    }
    b = (2 * datadim->nrefloci) - a;
    a += prior;
    b += prior;
    h = gsl_ran_beta(r, a, b); // NOTE NEEDED TO FLIP A AND B
    //h = gsl_ran_beta(r, b, a); // NOTE NEEDED TO FLIP A AND B
    //cerr << n << ": " << a << ", " << b << ", " << h << endl;
    gsl_histogram_increment(param->hihist[n].hist, h);
    gsl_matrix_set(param->hi, n, 0, h);
  }
}

/* set phi to hybrid index*/
void setPhiHRef(paramcont * param, datadimcont * datadim, int N){
  int i, n;
  double onePhi;
  
  for(i=0; i<N; i++){
    for(n=0; n<datadim->nind; n++){
      onePhi = gsl_matrix_get(param->hi, n, 0);
      gsl_matrix_set(param->refphi,i,n,onePhi);
    }
  }
}
