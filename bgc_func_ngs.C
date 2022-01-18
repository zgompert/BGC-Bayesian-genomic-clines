#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <float.h>

#include <math.h>

#include "bgc.h"

using namespace std;

/* -------------------------------------------------------------- */
/*                     functions for ngs data                     */
/* -------------------------------------------------------------- */

/* read in genotype likelihood data for the admixed and parental populations */
void getDatangs(string hybridfile, string p0file, string p1file, datacont * data,
		datadimcont * datadim){

  // read and store hybrid data
  readNgsFile(hybridfile,data->glh,datadim->nloci,datadim->nind);
  
  // read and store p0 data
  readNgsFile(p0file,data->glp0,datadim->nloci,datadim->nindP0);
  
  // read and store p1 data
  readNgsFile(p1file,data->glp1,datadim->nloci,datadim->nindP1);
}  

/* get key dimensions, nloci, ninds for gl format ngs data */
void getNgsDims(string hybridfile, string p0file, string p1file, datacont * data,
		datadimcont * datadim){

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
  datadim->nloci = atoi(element.c_str());  
  stream >> element; // number of inds
  datadim->nind = atoi(element.c_str());
  infile.close();

  // parent 0
  infile.open(p0file.c_str());
  if (!infile){
    cerr << "Cannot open file " << p0file << endl;
    exit(1);
  }
  // read line with data dimensions
  getline(infile, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of  loci
  stream >> element; // number of  individuals
  datadim->nindP0 = atoi(element.c_str());  
  infile.close();

  // parent 1
  infile.open(p1file.c_str());
  if (!infile){
    cerr << "Cannot open file " << p1file << endl;
    exit(1);
  }
  // read line with data dimensions
  getline(infile, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of  loci
  stream >> element; // number of  individuals
  datadim->nindP1 = atoi(element.c_str());  
  infile.close();
  
}

/* read gl format data from file */
void readNgsFile(string file, mat_container * gl, int L, int N){
  int i, j, a;
  double gensum = 0;
  double genotypes[GEN];
  string line, element;
  ifstream infile;
  istringstream stream;
  
  // open gl file
  infile.open(file.c_str());
  if (!infile){
    cerr << "Cannot open file " << file << endl;
    exit(1);
  }

  // read and burn dimensions header and list of individuals
  getline(infile, line);
  getline(infile, line);

  // read through and store data
  for(i=0; i<L; i++){ // loop over loci
    getline(infile, line); // data for one locus
    stream.str(line);
    stream.clear();
    stream >> element; // locus id, this is not retained
    for(j=0; j<N; j++){
      gensum = 0;
      for(a=0; a<GEN; a++){ // store genotype likelihoods
        stream >> element; // need to convert from phred scale: phred
        // = -10 log10(prob), prob = 10^(phred/-10)
        genotypes[a] = pow(10, (atof(element.c_str())/-10.0));
        gensum += genotypes[a];
      }
      for(a=0; a<GEN; a++){ // normalize genotype likelihoods
        genotypes[a] = genotypes[a]/gensum;
        if (genotypes[a] == 0) // set to DBL_MIN if 0
          genotypes[a] = DBL_MIN;
	gsl_matrix_set(gl[i].mat, j, a, genotypes[a]);
      }
    }
  }
}

/* function to update z, ancestry of hybrids, uses gibbs sampling similar to structure, but with ngs data */
/* we analytically sum across the unkown genotype */
void updateZngs(paramcont * param, datacont * data, datadimcont * datadim, auxcont * auxvar){
  int i, n, k, l;
  double p00 = 0, p01 = 0, p10 = 0, p11 = 0, ptot = 0;
  double pr[4]; // array of probabilities
  unsigned int sam[4]; // sampled ancestry for a locus
  int finite;
  double probtemp;

  // loop through loci
  for(i=0; i<datadim->nloci; i++){

    // loop through individuals
    for(n=0; n<datadim->nind; n++){
      p00 = 0;
      p01 = 0;
      p11 = 0;
      // first homozygous genotype
      k = 0; l = 0;
      // probability of data given the genotype, i.e., p(d | g) or gl
      probtemp = gsl_matrix_get(data->glh[i].mat, n, 0);
      // multiply by probability of genotype given ancestry
      // p(d | g) * p(g | z), note k = l
      p00 += probtemp * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l);
      p01 += probtemp * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
      p11 += probtemp * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);

      //cerr << i << "," << n << ": " << p00 << "," << p01 << "," << p11 << endl;
      
      // second homozygous genotype
      k = 1; l = 1;
      // probability of data given the genotype, i.e., p(d | g) or gl
      probtemp = gsl_matrix_get(data->glh[i].mat, n, 2);
      // multiply by probability of genotype given ancestry
      // p(d | g) * p(g | z), note k = l
      p00 += probtemp * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l);
      p01 += probtemp * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
      p11 += probtemp * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);

      //cerr << i << "," << n << ": " << p00 << "," << p01 << "," << p11 << endl;

      
      // heterozgyous genotype
      k = 0; l = 1;
      // probability of data given the genotype, i.e., p(d | g) or gl
      probtemp = gsl_matrix_get(data->glh[i].mat, n, 1);
      p00 += probtemp * 2 * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l);
      p01 += probtemp * (gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l) +
      			 gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l));
      p11 += probtemp * 2 * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);

      //cerr << i << "," << n << ": " << p00 << "," << p01 << "," << p11 << endl;
      
      // calculate p(z|phi) and mulitply by p(x|g) p(g|z,pi)
      p00 *= gsl_pow_2(1 - gsl_matrix_get(param->phi, i, n));
      p01 *= (1 - gsl_matrix_get(param->phi, i, n)) * gsl_matrix_get(param->phi, i, n);
      p11 *= gsl_pow_2(gsl_matrix_get(param->phi, i, n));

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
      if(sam[0] == 1){ // homo pop 0
	gsl_matrix_int_set(param->zCur[i].mat, n, 0, 0);
	gsl_matrix_int_set(param->zCur[i].mat, n, 1, 0);
      }
      else if(sam[1] == 1){ // het
	gsl_matrix_int_set(param->zCur[i].mat, n, 0, 0);
	gsl_matrix_int_set(param->zCur[i].mat, n, 1, 1);
      }
      else if(sam[2] == 1){ // het
	gsl_matrix_int_set(param->zCur[i].mat, n, 0, 1);
	gsl_matrix_int_set(param->zCur[i].mat, n, 1, 0);
      }
      else{ // homo pop 1
	gsl_matrix_int_set(param->zCur[i].mat, n, 0, 1);
	gsl_matrix_int_set(param->zCur[i].mat, n, 1, 1);
      }
    }
  }
}

/* initialize ancestry, this only works for ngs */
void initZngs(paramcont * param, datacont * data, datadimcont * datadim, auxcont * auxvar){
  int i, n, k, l;
  double p00 = 0, p01 = 0, p10 = 0, p11 = 0, ptot = 0;
  double pr[4]; // array of probabilities
  unsigned int sam[4]; // sampled ancestry for a locus
  double probtemp;
  int finite;

  // loop through loci
  for(i=0; i<datadim->nloci; i++){

    // loop through individuals
    for(n=0; n<datadim->nind; n++){
      p00 = 0;
      p01 = 0;
      p11 = 0;
      // first homozygous genotype
      k = 0; l = 0;
      // probability of data given the genotype, i.e., p(d | g) or gl
      probtemp = gsl_matrix_get(data->glh[i].mat, n, 0);
      // multiply by probability of genotype given ancestry
      // p(d | g) * p(g | z), note k = l
      p00 += probtemp * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l);
      p01 += probtemp * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
      p11 += probtemp * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);

      // second homozygous genotype
      k = 1; l = 1;
      // probability of data given the genotype, i.e., p(d | g) or gl
      probtemp = gsl_matrix_get(data->glh[i].mat, n, 2);
      // multiply by probability of genotype given ancestry
      // p(d | g) * p(g | z), note k = l
      p00 += probtemp * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l);
      p01 += probtemp * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
      p11 += probtemp * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);

      // heterozgyous genotype
      k = 0; l = 1;
      // probability of data given the genotype, i.e., p(d | g) or gl
      probtemp = gsl_matrix_get(data->glh[i].mat, n, 1);
      // multiply by probability of genotype given ancestry
      // p(d | g) * p(g | z), note k = l
      p00 += probtemp * 2 * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l);
      p01 += probtemp * (gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l) +
			 gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l));
      p11 += probtemp * 2 * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
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
      // sample ancestry pair based on probabilities
      gsl_ran_multinomial(r, 4, 1, pr, sam);
      if(sam[0] == 1){ // homo pop 0
	gsl_matrix_int_set(param->zCur[i].mat, n, 0, 0);
	gsl_matrix_int_set(param->zCur[i].mat, n, 1, 0);
      }
      else if(sam[1] == 1){ // het
	gsl_matrix_int_set(param->zCur[i].mat, n, 0, 0);
	gsl_matrix_int_set(param->zCur[i].mat, n, 1, 1);
      }
      else if(sam[2] == 1){ // het
	gsl_matrix_int_set(param->zCur[i].mat, n, 0, 1);
	gsl_matrix_int_set(param->zCur[i].mat, n, 1, 0);
      }
      else{ // homo pop 1
	gsl_matrix_int_set(param->zCur[i].mat, n, 0, 1);
	gsl_matrix_int_set(param->zCur[i].mat, n, 1, 1);
      }
    }
  }
}

// sample genotype based on gl and add to count
void initPcnt(gsl_matrix_int * pcnt, mat_container * gl, int nloci, int nind){

  int i, k, n;
  double prob[GEN];
  unsigned int geno[GEN];
  
  
  // set pcnt to zero
  gsl_matrix_int_set_zero(pcnt);

  for(i=0; i<nloci; i++){
    for(n=0; n<nind; n++){
      for(k=0; k<GEN; k++){
	prob[k] = gsl_matrix_get(gl[i].mat, n, k);
      }
      gsl_ran_multinomial(r, GEN, 1, prob, geno);
      for(k=0; k<GEN; k++){
	if(geno[k] == 1){
	  gsl_matrix_int_set(pcnt, i, 0, k);
	  gsl_matrix_int_set(pcnt, i, 1, 2-k); // two columns, one for each allele
	}
      }
    }
  }
}

/* gibbs sample update of parental population counts */
void updatePcnt(datacont * data, auxcont * auxvar, datadimcont * datadim, gsl_matrix_int * pcnt, 
		mat_container * gl, int nind, gsl_matrix * pi){
  int i, k, n, a, b;
  double prob[GEN];
  unsigned int geno[GEN];

  // set pcnt to zero
  gsl_matrix_int_set_zero(pcnt);

  for(i=0; i<datadim->nloci; i++){
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
    //cerr << gsl_matrix_int_get(pcnt, i, 0) << ", ";
  }
  //cerr << endl;
}


/* calculate Ln Likelihood of the data from the admixed populations = p(x | z, pi) = \prod_i,n sum_g p(x | g) p(g | pi, z) */
// double calcLnLngs(paramcont * param, datacont * data, datadimcont * datadim, auxcont * auxvar){

//   int i, n, k, l, x, nk;
//   int acopy = 0;
//   double probsum, probtemp, prob = 0;
//   double err;
  
//   // loop through loci
//   for(i=0; i<datadim->nloci; i++){
//     data->error = gsl_vector_get(data->errVec, i);
//     nk = gsl_vector_int_get(datadim->nallele, i);
//     err = (double) data->error / (nk - 1);
//     // loop through individuals
//     for(n=0; n<datadim->nind; n++){
//       probsum = 0;
//       acopy = 0;
//       for(k=0; k<nk; k++){
// 	gsl_vector_uint_set(auxvar->vecAllele_int, k,
// 			    gsl_matrix_int_get(data->phcnt[i].mat, n, k));
// 	// the individual has this allele
// 	if(gsl_matrix_int_get(data->phcnt[i].mat, n, k) > 0){
// 	  acopy++;
// 	}
//       }
//       // if acopy == 0 no data for this locus, do not include in likelihood calculation
//       if (acopy > 0){
// 	for(k=0; k<nk; k++){
// 	  for(l=k; l<nk; l++){
// 	    // this is a homozygous genotype
// 	    if (k == l){
// 	      // set multinomial probs. for homozygous 
// 	      for (x=0; x<nk; x++){
// 		if (k == x){
// 		  gsl_vector_set(auxvar->vecAllelen1, x, (1.0 - data->error));
// 		}
// 		else {
// 		  gsl_vector_set(auxvar->vecAllelen1, x, err);
// 		}
// 	      } 
// 	      // probability of data given the genotype, i.e., p(x | g)
// 	      probtemp = gsl_ran_multinomial_pdf(nk, auxvar->vecAllelen1->data, auxvar->vecAllele_int->data);
// 	      // probability of genotype | ancestry
// 	      // ancestry only pop. 0
// 	      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
// 		  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
// 		probtemp = probtemp * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l);
// 	      }
// 	      // ancestry only pop. 1	      
// 	      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
// 		       gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
// 		probtemp = probtemp * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
// 	      }
// 	      // ancestry from each population 
// 	      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) != gsl_matrix_int_get(param->zCur[i].mat, n, 1)){	      
// 		probtemp = probtemp * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
// 	      }
// 	      if (probtemp == 0){
// 		probtemp = DBL_MIN;
// 	      }
// 	      probsum += probtemp;
// 	    }
// 	    // this is a heterozygous genotype
// 	    else { 
// 	      // set multinomial probs. for homozygous 
// 	      for (x=0; x<nk; x++){
// 		if (x == k || x == l){
// 		  gsl_vector_set(auxvar->vecAllelen1, x, (0.5 - data->error + 0.5 * err));
// 		}
// 		else {
// 		  gsl_vector_set(auxvar->vecAllelen1, x, err);
// 		}
// 	      } 
// 	      // probability of data given the genotype, i.e., p(x | g)
// 	      probtemp = gsl_ran_multinomial_pdf(nk, auxvar->vecAllelen1->data, auxvar->vecAllele_int->data);
// 	      // probability of genotype | ancestry
// 	      // ancestry only pop. 0
// 	      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
// 		  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
// 		probtemp = probtemp * 2 * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l);
// 	      }
// 	      // ancestry only pop. 1	      
// 	      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
// 		       gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
// 		probtemp = probtemp * 2 * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
// 	      }
// 	      // ancestry from each population 
// 	      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) != gsl_matrix_int_get(param->zCur[i].mat, n, 1)){	      
// 		probtemp =  probtemp * (gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l) +
// 					gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l));
// 	      }
// 	      if (probtemp == 0){
// 		probtemp = DBL_MIN;
// 	      }
// 	      probsum += probtemp;
// 	    }
// 	  }
// 	}
// 	prob += log(probsum);
//       }
//     }
//   }
//   return(prob);
// }
