//
// File: RGBColor.h
// Created by: Julien Dutheil
// Created on: Thu Mar 16 2006
//

/*
Copyright or � or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide utilitary
classes. This file belongs to the Bio++ Project.

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

#ifndef _RGBCOLOR_H_
#define _RGBCOLOR_H_

// From the STL:
#include <cstdlib>
#include <cmath>
using namespace std;

#include "TextTools.h"
#include "Exceptions.h"
#include "Clonable.h"

namespace bpp
{

/**
 * @brief Describe a color according to its red, green and blue componants.
 */
class RGBColor:
  public virtual Clonable
{
  protected:
    unsigned int _red;
    unsigned int _green;
    unsigned int _blue;

  public:
    RGBColor(unsigned int red, unsigned int green, unsigned int blue): _red(red), _green(green), _blue(blue) {} 
    RGBColor(): _red(0), _green(0), _blue(0) {} 
    virtual ~RGBColor() {}

#ifdef NO_VIRTUAL_COV
    Clonable*
#else
    RGBColor*
#endif
    clone() const { return new RGBColor(*this); }

  public:
    bool operator==(const RGBColor & color) const
    {
      return _red == color._red && _green == color._green && _blue == color._blue;
    }

    /**
     * @brief Comparison operator (for sorting purposes).
     *
     * The hexadecimal string representation is used for comparison.
     *
     * @param color The color to compare with.
     */
    bool operator<(const RGBColor & color) const
    {
      return toHex() < color.toHex();
    }

    /**
     * @brief Get the HTML-like, hexadecimal description of this color.
     */
    string toHex() const
    {
      string hex = "#";
      hex += decToHex(_red);
      hex += decToHex(_green);
      hex += decToHex(_blue);
      return hex;
    }

    /**
     * @brief Access to each color componant: 0=red, 1=green, 2=blue.
     */
    const unsigned int & operator[](unsigned int i) const
    {
      if(i == 0) return _red;
      if(i == 1) return _green;
      if(i == 2) return _blue;
      throw Exception("Invalid color index");
    }

    /**
     * @brief Access to each color componant: 0=red, 1=green, 2=blue.
     */
    unsigned int & operator[](unsigned int i)
    {
      if(i == 0) return _red;
      if(i == 1) return _green;
      if(i == 2) return _blue;
      throw Exception("Invalid color index");
    }

    /**
     * @brief Get a string description of the color, e.g. [R255,G0,B255].
     */
    string toString() const
    {
      return "[R" + TextTools::toString(_red) + ",G" + TextTools::toString(_green) + ",B" + TextTools::toString(_blue) + "]"; 
    }

  protected:
    static string decToHex(unsigned int dec)
    {
      string hexa = "0123456789ABCDEF";
      string hex = "";
      while (dec > 15)
      {
        unsigned int tmp = dec - (int)floor((double)dec/16.)*16;
        hex = hexa[tmp] + hex;
        dec=(int)floor((double)dec/16.);
      }
      hex = hexa[dec] + hex;
      if(hex.size() == 1) hex = "0" + hex;
      return hex;
    }

};

} // end of namespace bpp;

#endif //_RGBCOLOR_H_

