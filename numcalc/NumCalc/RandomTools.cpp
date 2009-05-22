//
// File RandomTools.cpp
// Author : Julien Dutheil
// Last modification : Friday Septembre 24 2004
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "RandomTools.h"
#include "VectorTools.h"
#include "Uniform01K.h"

using namespace bpp;

RandomFactory * RandomTools::DEFAULT_GENERATOR = new Uniform01K(time(NULL));

// Initiate random seed :
//RandomTools::RandInt RandomTools::r = time(NULL) ;

void RandomTools::setSeed(long seed)
{
  DEFAULT_GENERATOR->setSeed(seed);
}

// Method to get a double random value (between 0 and specified range)
// Note : the number you get is between 0 and entry not including entry !
double RandomTools::giveRandomNumberBetweenZeroAndEntry(double entry, const RandomFactory * generator)
{
  //double tm = r.drawFloatNumber();
  double tm = generator->drawNumber();
  return (tm * entry);
}

// Method to get a boolean random value
bool RandomTools::flipCoin(const RandomFactory * generator)
{
  return ((RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0, generator) - 0.5) > 0);
}

// Method to get a integer random value (between 0 and specified range)
// Note : the number you get is between 0 and entry not including entry !
int RandomTools::giveIntRandomNumberBetweenZeroAndEntry(int entry, const RandomFactory * generator)
{
  return (int)(giveRandomNumberBetweenZeroAndEntry(entry, generator));
}

double RandomTools::randGaussian(double mean, double variance, const RandomFactory * generator)
{
  static int N = 100;
  static double X;
  X=0.0-N/2; /* set mean to 0 */
  for (int ri = 0; ri < N; ri++)
  {
    //    X += 1.0*rand()/RAND_MAX;
    X += giveRandomNumberBetweenZeroAndEntry(1, generator);
  }
  
  /* for uniform randoms in [0,1], mu = 0.5 and var = 1/12 */
  /* adjust X so mu = 0 and var = 1 */
    
  //  X = X * sqrt(12 / N);       /* adjust variance to 1 */
  //  cout <<X * sqrt(variance*12.0/N) + mean<<" ";
  double g = X * sqrt(variance*12.0/N) + mean;
  return (g);
}

double RandomTools::randGamma(double dblAlpha, const RandomFactory * generator)
{
  assert(dblAlpha > 0.0);
  if( dblAlpha < 1.0 ) return RandomTools::DblGammaLessThanOne(dblAlpha, generator);
  else if( dblAlpha > 1.0 ) return RandomTools::DblGammaGreaterThanOne(dblAlpha, generator);
  return -log(RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0, generator));
}  

double RandomTools::randGamma(double alpha, double beta, const RandomFactory * generator)
{
  double x= RandomTools::randGamma(alpha, generator) / beta;
  return x;
}

double RandomTools::randExponential(double mean, const RandomFactory * generator)
{
  return - mean * log(RandomTools::giveRandomNumberBetweenZeroAndEntry(1, generator));
}

vector<unsigned int> RandomTools::randMultinomial(unsigned int n, const vector<double>& probs)
{
  double s = VectorTools::sum(probs);
  double r;
  double cumprob;
  vector<unsigned int> sample(n);
  for(unsigned int i = 0; i < n; i++)
  {
    r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
    cumprob = 0;
    bool test = true;
    for(unsigned int j = 0; test & (j < probs.size()); j++)
    {
      cumprob += probs[j] / s;
      if(r <= cumprob)
      {
        sample[i] = j;
        test = false;
      }
    }
    // This test should never be true if probs sum to one:
    if(test) sample[i] = probs.size();
  }
  return sample;
}

//------------------------------------------------------------------------------

  
double RandomTools::DblGammaGreaterThanOne(double dblAlpha, const RandomFactory * generator)
{
  // Code adopted from David Heckerman
  //-----------------------------------------------------------
  //  DblGammaGreaterThanOne(dblAlpha)
  //
  //  routine to generate a gamma random variable with unit scale and
  //      alpha > 1
  //  reference: Ripley, Stochastic Simulation, p.90 
  //  Chang and Feast, Appl.Stat. (28) p.290
  //-----------------------------------------------------------
  double rgdbl[6];
    
  rgdbl[1] = dblAlpha - 1.0;
  rgdbl[2] = (dblAlpha - (1.0 / (6.0 * dblAlpha))) / rgdbl[1];
  rgdbl[3] = 2.0 / rgdbl[1];
  rgdbl[4] = rgdbl[3] + 2.0;
  rgdbl[5] = 1.0 / sqrt(dblAlpha);
    
  for (;;)
  {
    double dblRand1;
    double dblRand2;
    do
    {
      dblRand1 = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0, generator);
      dblRand2 = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0, generator);
      if (dblAlpha > 2.5) dblRand1 = dblRand2 + rgdbl[5] * (1.0 - 1.86 * dblRand1);
    }
    while (!(0.0 < dblRand1 && dblRand1 < 1.0));
  
    double dblTemp = rgdbl[2] * dblRand2 / dblRand1;
  
    if (rgdbl[3] * dblRand1 + dblTemp + 1.0 / dblTemp <= rgdbl[4] ||
        rgdbl[3] * log(dblRand1) + dblTemp - log(dblTemp) < 1.0)
    {
      return dblTemp * rgdbl[1];
    }
  }
  assert(false);
  return 0.0;
}

double RandomTools::DblGammaLessThanOne(double dblAlpha, const RandomFactory * generator) {
  //routine to generate a gamma random variable with 
  //unit scale and alpha < 1
  //reference: Ripley, Stochastic Simulation, p.88 
  double dblTemp;
  const double dblexp = exp(1.0);
  for (;;)
  {
    double dblRand0 = giveRandomNumberBetweenZeroAndEntry(1.0, generator);
    double dblRand1 = giveRandomNumberBetweenZeroAndEntry(1.0, generator);
    if (dblRand0 <= (dblexp / (dblAlpha + dblexp)))
    {
      dblTemp = pow(((dblAlpha + dblexp) * dblRand0) /
      dblexp, 1.0 / dblAlpha);
      if (dblRand1 <= exp(-1.0 * dblTemp)) return dblTemp;
    }
    else
    {
      dblTemp = -1.0 * log((dblAlpha + dblexp) * (1.0 - dblRand0) / (dblAlpha * dblexp)); 
      if (dblRand1 <= pow(dblTemp,dblAlpha - 1.0)) return dblTemp;
    }
  }
  assert(false);
  return 0.0;
}

/******************************************************************************/

//From Yang's PAML package:

/******************************************************************************/

double RandomTools::qNorm(double prob)
{
   double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   double y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) return (-9999);

   y = sqrt (log(1/(p1*p1)));   
   z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   return (p<0.5 ? -z : z);
}


double RandomTools::lnGamma (double alpha)
{
  double x=alpha, f=0, z;

  if (x<7)
  {
    f=1;  z=x-1;
    while (++z<7)  f*=z;
    x=z;   f=-log(f);
  }
  z = 1/(x*x);
  return  f + (x-0.5)*log(x) - x + .918938533204673 
    + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
         +.083333333333333)/x;  
}



double RandomTools::incompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
  int i;
  double p=alpha, g=ln_gamma_alpha;
  double accurate=1e-8, overflow=1e30;
  double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0;
  vector<double> pn(6);

  if (x==0) return (0);
  if (x<0 || p<=0) return (-1);

  factor=exp(p*log(x)-x-g);   
  if (x>1 && x>=p) goto l30;
  /* (1) series expansion */
  gin=1;  term=1;  rn=p;
l20:
  rn++;
  term*=x/rn;   gin+=term;

  if (term > accurate) goto l20;
  gin*=factor/p;
  goto l50;
l30:
  /* (2) continued fraction */
  a=1-p;   b=a+x+1;  term=0;
  pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
  gin=pn[2]/pn[3];
l32:
  a++;  b+=2;  term++;   an=a*term;
  for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
  if (pn[5] == 0) goto l35;
  rn=pn[4]/pn[5];   dif=fabs(gin-rn);
  if (dif>accurate) goto l34;
  if (dif<=accurate*rn) goto l42;
l34:
  gin=rn;
l35:
  for (i=0; i<4; i++) pn[i]=pn[i+2];
  if (fabs(pn[4]) < overflow) goto l32;
  for (i=0; i<4; i++) pn[i]/=overflow;
  goto l32;
l42:
  gin=1-factor*gin;

l50:

  return (gin);
}



double RandomTools::qChisq(double prob, double v)
{
  double e=.5e-6, aa=.6931471805, p=prob, g;
  double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

  if (p<.000002 || p>.999998 || v<=0) return (-1);

  g = lnGamma (v/2);
  xx=v/2;   c=xx-1;
  if (v >= -1.24*log(p)) goto l1;

  ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
  if (ch-e<0) return (ch);
  goto l4;
l1:
  if (v>.32) goto l3;
  ch=0.4;   a=log(1-p);
l2:
  q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
  t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
  ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
  if (fabs(q/ch-1)-.01 <= 0) goto l4;
  else                       goto l2;
  
l3: 
  x=qNorm (p);
  p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
  if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
  q=ch;   p1=.5*ch;
  if ((t=incompleteGamma (p1, xx, g))<0)
  {
    printf ("\nerr IncompleteGamma");
    return (-1);
  }
  p2=p-t;
  t=p2*exp(xx*aa+g+p1-c*log(ch));   
  b=t/ch;  a=0.5*t-b*c;

  s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
  s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
  s3=(210+a*(462+a*(707+932*a)))/2520;
  s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
  s5=(84+264*a+c*(175+606*a))/2520;
  s6=(120+c*(346+127*c))/5040;
  ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
  if (fabs(q/ch-1) > e) goto l4;

  return (ch);
}

//------------------------------------------------------------------------------

