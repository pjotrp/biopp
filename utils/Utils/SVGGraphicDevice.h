//
// File: SVGGraphicDevice.h
// Created by: Julien Dutheil
// Created on: Mon Mar 10 2008
//

/*
Copyright or © or Copr. CNRS, (November 16, 2006)

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

#ifndef _SVGGRAPHICDEVICE_H_
#define _SVGGRAPHICDEVICE_H_

#include "GraphicDevice.h"
#include "ColorTools.h"

// From the STL:
#include <map>
using namespace std;

namespace bpp
{

/**
 * @brief SVG plotting format.
 */
class SVGGraphicDevice:
  public virtual GraphicDevice
{
  protected:
    ostream & _out;
    RGBColor _fgColor;
    RGBColor _bgColor;
    Font _font; 
    unsigned int _pointSize;
    short _lineType;
    int _currentLayer;
    map<int, vector<string>, greater<int> > _layers; //Layer display as in xfig
    bool _inkscapeEnabled;

  public:
    SVGGraphicDevice(ostream & out, bool inkscapeEnabled = false):
      _out(out),
      _fgColor(ColorTools::BLACK),
      _bgColor(ColorTools::WHITE),
      _font(),
      _pointSize(1),
      _inkscapeEnabled(inkscapeEnabled)
    {}

    virtual ~SVGGraphicDevice() {}

  public:
    void beginDocument();
    void endDocument();

    void setCurrentForegroundColor(const RGBColor & color);
    void setCurrentBackgroundColor(const RGBColor & color);
    void setCurrentFont(const Font & font);
    void setCurrentPointSize(unsigned int size) { _pointSize = size; }
    void setCurrentLineType(short type) throw (Exception)
    { 
      if(type == LINE_SOLID) _lineType = type;
      else if(type == LINE_DASHED) _lineType = type;
      else if(type == LINE_DOTTED) _lineType = type;
      else throw Exception("SVGGraphicDevice::setCurrentLineType. Unknown line type: " + TextTools::toString(type));
    }
    void setCurrentLayer(int layerIndex) { _currentLayer = layerIndex; }
    void drawLine(double x1, double y1, double x2, double y2);
    void drawRect(double x, double y, double width, double height, short fill = FILL_EMPTY);
    void drawCircle(double x, double y, double radius, short fill = FILL_EMPTY);
    void drawText(double x, double y, const string & text, short hpos = TEXT_HORIZONTAL_LEFT, short vpos = TEXT_VERTICAL_BOTTOM, double angle = 0);
    void comment(const string & comment)
    {
      _layers[_currentLayer].push_back("<!-- " + comment + " -->");
    }

  public:
    static string colorToText(const RGBColor & color)
    {
      return "rgb(" + TextTools::toString(color[0]) + "," + TextTools::toString(color[1]) + "," + TextTools::toString(color[2]) + ")";
    }

};

} // end of namespace bpp.

#endif //_SVGGRAPHICDEVICE_H_


