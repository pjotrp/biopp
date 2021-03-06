//
// File: RColorSet.cpp
// Created by: Julien Dutheil
// Created on: Mon Apr 14 2008
//

/*
Copyright or © or Copr. CNRS, (November 17, 2008)

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

#include "RColorSet.h"

using namespace bpp;

RColorSet::RColorSet()
{
  _colors["white"] = RGBColor(255, 255, 255);
  _colors["aliceblue"] = RGBColor(240, 248, 255);
  _colors["antiquewhite"] = RGBColor(250, 235, 215);
  _colors["antiquewhite1"] = RGBColor(255, 239, 219);
  _colors["antiquewhite2"] = RGBColor(238, 223, 204);
  _colors["antiquewhite3"] = RGBColor(205, 192, 176);
  _colors["antiquewhite4"] = RGBColor(139, 131, 120);
  _colors["aquamarine"] = RGBColor(127, 255, 212);
  _colors["aquamarine1"] = RGBColor(127, 255, 212);
  _colors["aquamarine2"] = RGBColor(118, 238, 198);
  _colors["aquamarine3"] = RGBColor(102, 205, 170);
  _colors["aquamarine4"] = RGBColor(69, 139, 116);
  _colors["azure"] = RGBColor(240, 255, 255);
  _colors["azure1"] = RGBColor(240, 255, 255);
  _colors["azure2"] = RGBColor(224, 238, 238);
  _colors["azure3"] = RGBColor(193, 205, 205);
  _colors["azure4"] = RGBColor(131, 139, 139);
  _colors["beige"] = RGBColor(245, 245, 220);
  _colors["bisque"] = RGBColor(255, 228, 196);
  _colors["bisque1"] = RGBColor(255, 228, 196);
  _colors["bisque2"] = RGBColor(238, 213, 183);
  _colors["bisque3"] = RGBColor(205, 183, 158);
  _colors["bisque4"] = RGBColor(139, 125, 107);
  _colors["black"] = RGBColor(0, 0, 0);
  _colors["blanchedalmond"] = RGBColor(255, 235, 205);
  _colors["blue"] = RGBColor(0, 0, 255);
  _colors["blue1"] = RGBColor(0, 0, 255);
  _colors["blue2"] = RGBColor(0, 0, 238);
  _colors["blue3"] = RGBColor(0, 0, 205);
  _colors["blue4"] = RGBColor(0, 0, 139);
  _colors["blueviolet"] = RGBColor(138, 43, 226);
  _colors["brown"] = RGBColor(165, 42, 42);
  _colors["brown1"] = RGBColor(255, 64, 64);
  _colors["brown2"] = RGBColor(238, 59, 59);
  _colors["brown3"] = RGBColor(205, 51, 51);
  _colors["brown4"] = RGBColor(139, 35, 35);
  _colors["burlywood"] = RGBColor(222, 184, 135);
  _colors["burlywood1"] = RGBColor(255, 211, 155);
  _colors["burlywood2"] = RGBColor(238, 197, 145);
  _colors["burlywood3"] = RGBColor(205, 170, 125);
  _colors["burlywood4"] = RGBColor(139, 115, 85);
  _colors["cadetblue"] = RGBColor(95, 158, 160);
  _colors["cadetblue1"] = RGBColor(152, 245, 255);
  _colors["cadetblue2"] = RGBColor(142, 229, 238);
  _colors["cadetblue3"] = RGBColor(122, 197, 205);
  _colors["cadetblue4"] = RGBColor(83, 134, 139);
  _colors["chartreuse"] = RGBColor(127, 255, 0);
  _colors["chartreuse1"] = RGBColor(127, 255, 0);
  _colors["chartreuse2"] = RGBColor(118, 238, 0);
  _colors["chartreuse3"] = RGBColor(102, 205, 0);
  _colors["chartreuse4"] = RGBColor(69, 139, 0);
  _colors["chocolate"] = RGBColor(210, 105, 30);
  _colors["chocolate1"] = RGBColor(255, 127, 36);
  _colors["chocolate2"] = RGBColor(238, 118, 33);
  _colors["chocolate3"] = RGBColor(205, 102, 29);
  _colors["chocolate4"] = RGBColor(139, 69, 19);
  _colors["coral"] = RGBColor(255, 127, 80);
  _colors["coral1"] = RGBColor(255, 114, 86);
  _colors["coral2"] = RGBColor(238, 106, 80);
  _colors["coral3"] = RGBColor(205, 91, 69);
  _colors["coral4"] = RGBColor(139, 62, 47);
  _colors["cornflowerblue"] = RGBColor(100, 149, 237);
  _colors["cornsilk"] = RGBColor(255, 248, 220);
  _colors["cornsilk1"] = RGBColor(255, 248, 220);
  _colors["cornsilk2"] = RGBColor(238, 232, 205);
  _colors["cornsilk3"] = RGBColor(205, 200, 177);
  _colors["cornsilk4"] = RGBColor(139, 136, 120);
  _colors["cyan"] = RGBColor(0, 255, 255);
  _colors["cyan1"] = RGBColor(0, 255, 255);
  _colors["cyan2"] = RGBColor(0, 238, 238);
  _colors["cyan3"] = RGBColor(0, 205, 205);
  _colors["cyan4"] = RGBColor(0, 139, 139);
  _colors["darkblue"] = RGBColor(0, 0, 139);
  _colors["darkcyan"] = RGBColor(0, 139, 139);
  _colors["darkgoldenrod"] = RGBColor(184, 134, 11);
  _colors["darkgoldenrod1"] = RGBColor(255, 185, 15);
  _colors["darkgoldenrod2"] = RGBColor(238, 173, 14);
  _colors["darkgoldenrod3"] = RGBColor(205, 149, 12);
  _colors["darkgoldenrod4"] = RGBColor(139, 101, 8);
  _colors["darkgray"] = RGBColor(169, 169, 169);
  _colors["darkgreen"] = RGBColor(0, 100, 0);
  _colors["darkgrey"] = RGBColor(169, 169, 169);
  _colors["darkkhaki"] = RGBColor(189, 183, 107);
  _colors["darkmagenta"] = RGBColor(139, 0, 139);
  _colors["darkolivegreen"] = RGBColor(85, 107, 47);
  _colors["darkolivegreen1"] = RGBColor(202, 255, 112);
  _colors["darkolivegreen2"] = RGBColor(188, 238, 104);
  _colors["darkolivegreen3"] = RGBColor(162, 205, 90);
  _colors["darkolivegreen4"] = RGBColor(110, 139, 61);
  _colors["darkorange"] = RGBColor(255, 140, 0);
  _colors["darkorange1"] = RGBColor(255, 127, 0);
  _colors["darkorange2"] = RGBColor(238, 118, 0);
  _colors["darkorange3"] = RGBColor(205, 102, 0);
  _colors["darkorange4"] = RGBColor(139, 69, 0);
  _colors["darkorchid"] = RGBColor(153, 50, 204);
  _colors["darkorchid1"] = RGBColor(191, 62, 255);
  _colors["darkorchid2"] = RGBColor(178, 58, 238);
  _colors["darkorchid3"] = RGBColor(154, 50, 205);
  _colors["darkorchid4"] = RGBColor(104, 34, 139);
  _colors["darkred"] = RGBColor(139, 0, 0);
  _colors["darksalmon"] = RGBColor(233, 150, 122);
  _colors["darkseagreen"] = RGBColor(143, 188, 143);
  _colors["darkseagreen1"] = RGBColor(193, 255, 193);
  _colors["darkseagreen2"] = RGBColor(180, 238, 180);
  _colors["darkseagreen3"] = RGBColor(155, 205, 155);
  _colors["darkseagreen4"] = RGBColor(105, 139, 105);
  _colors["darkslateblue"] = RGBColor(72, 61, 139);
  _colors["darkslategray"] = RGBColor(47, 79, 79);
  _colors["darkslategray1"] = RGBColor(151, 255, 255);
  _colors["darkslategray2"] = RGBColor(141, 238, 238);
  _colors["darkslategray3"] = RGBColor(121, 205, 205);
  _colors["darkslategray4"] = RGBColor(82, 139, 139);
  _colors["darkslategrey"] = RGBColor(47, 79, 79);
  _colors["darkturquoise"] = RGBColor(0, 206, 209);
  _colors["darkviolet"] = RGBColor(148, 0, 211);
  _colors["deeppink"] = RGBColor(255, 20, 147);
  _colors["deeppink1"] = RGBColor(255, 20, 147);
  _colors["deeppink2"] = RGBColor(238, 18, 137);
  _colors["deeppink3"] = RGBColor(205, 16, 118);
  _colors["deeppink4"] = RGBColor(139, 10, 80);
  _colors["deepskyblue"] = RGBColor(0, 191, 255);
  _colors["deepskyblue1"] = RGBColor(0, 191, 255);
  _colors["deepskyblue2"] = RGBColor(0, 178, 238);
  _colors["deepskyblue3"] = RGBColor(0, 154, 205);
  _colors["deepskyblue4"] = RGBColor(0, 104, 139);
  _colors["dimgray"] = RGBColor(105, 105, 105);
  _colors["dimgrey"] = RGBColor(105, 105, 105);
  _colors["dodgerblue"] = RGBColor(30, 144, 255);
  _colors["dodgerblue1"] = RGBColor(30, 144, 255);
  _colors["dodgerblue2"] = RGBColor(28, 134, 238);
  _colors["dodgerblue3"] = RGBColor(24, 116, 205);
  _colors["dodgerblue4"] = RGBColor(16, 78, 139);
  _colors["firebrick"] = RGBColor(178, 34, 34);
  _colors["firebrick1"] = RGBColor(255, 48, 48);
  _colors["firebrick2"] = RGBColor(238, 44, 44);
  _colors["firebrick3"] = RGBColor(205, 38, 38);
  _colors["firebrick4"] = RGBColor(139, 26, 26);
  _colors["floralwhite"] = RGBColor(255, 250, 240);
  _colors["forestgreen"] = RGBColor(34, 139, 34);
  _colors["gainsboro"] = RGBColor(220, 220, 220);
  _colors["ghostwhite"] = RGBColor(248, 248, 255);
  _colors["gold"] = RGBColor(255, 215, 0);
  _colors["gold1"] = RGBColor(255, 215, 0);
  _colors["gold2"] = RGBColor(238, 201, 0);
  _colors["gold3"] = RGBColor(205, 173, 0);
  _colors["gold4"] = RGBColor(139, 117, 0);
  _colors["goldenrod"] = RGBColor(218, 165, 32);
  _colors["goldenrod1"] = RGBColor(255, 193, 37);
  _colors["goldenrod2"] = RGBColor(238, 180, 34);
  _colors["goldenrod3"] = RGBColor(205, 155, 29);
  _colors["goldenrod4"] = RGBColor(139, 105, 20);
  _colors["gray"] = RGBColor(190, 190, 190);
  _colors["gray0"] = RGBColor(0, 0, 0);
  _colors["gray1"] = RGBColor(3, 3, 3);
  _colors["gray2"] = RGBColor(5, 5, 5);
  _colors["gray3"] = RGBColor(8, 8, 8);
  _colors["gray4"] = RGBColor(10, 10, 10);
  _colors["gray5"] = RGBColor(13, 13, 13);
  _colors["gray6"] = RGBColor(15, 15, 15);
  _colors["gray7"] = RGBColor(18, 18, 18);
  _colors["gray8"] = RGBColor(20, 20, 20);
  _colors["gray9"] = RGBColor(23, 23, 23);
  _colors["gray10"] = RGBColor(26, 26, 26);
  _colors["gray11"] = RGBColor(28, 28, 28);
  _colors["gray12"] = RGBColor(31, 31, 31);
  _colors["gray13"] = RGBColor(33, 33, 33);
  _colors["gray14"] = RGBColor(36, 36, 36);
  _colors["gray15"] = RGBColor(38, 38, 38);
  _colors["gray16"] = RGBColor(41, 41, 41);
  _colors["gray17"] = RGBColor(43, 43, 43);
  _colors["gray18"] = RGBColor(46, 46, 46);
  _colors["gray19"] = RGBColor(48, 48, 48);
  _colors["gray20"] = RGBColor(51, 51, 51);
  _colors["gray21"] = RGBColor(54, 54, 54);
  _colors["gray22"] = RGBColor(56, 56, 56);
  _colors["gray23"] = RGBColor(59, 59, 59);
  _colors["gray24"] = RGBColor(61, 61, 61);
  _colors["gray25"] = RGBColor(64, 64, 64);
  _colors["gray26"] = RGBColor(66, 66, 66);
  _colors["gray27"] = RGBColor(69, 69, 69);
  _colors["gray28"] = RGBColor(71, 71, 71);
  _colors["gray29"] = RGBColor(74, 74, 74);
  _colors["gray30"] = RGBColor(77, 77, 77);
  _colors["gray31"] = RGBColor(79, 79, 79);
  _colors["gray32"] = RGBColor(82, 82, 82);
  _colors["gray33"] = RGBColor(84, 84, 84);
  _colors["gray34"] = RGBColor(87, 87, 87);
  _colors["gray35"] = RGBColor(89, 89, 89);
  _colors["gray36"] = RGBColor(92, 92, 92);
  _colors["gray37"] = RGBColor(94, 94, 94);
  _colors["gray38"] = RGBColor(97, 97, 97);
  _colors["gray39"] = RGBColor(99, 99, 99);
  _colors["gray40"] = RGBColor(102, 102, 102);
  _colors["gray41"] = RGBColor(105, 105, 105);
  _colors["gray42"] = RGBColor(107, 107, 107);
  _colors["gray43"] = RGBColor(110, 110, 110);
  _colors["gray44"] = RGBColor(112, 112, 112);
  _colors["gray45"] = RGBColor(115, 115, 115);
  _colors["gray46"] = RGBColor(117, 117, 117);
  _colors["gray47"] = RGBColor(120, 120, 120);
  _colors["gray48"] = RGBColor(122, 122, 122);
  _colors["gray49"] = RGBColor(125, 125, 125);
  _colors["gray50"] = RGBColor(127, 127, 127);
  _colors["gray51"] = RGBColor(130, 130, 130);
  _colors["gray52"] = RGBColor(133, 133, 133);
  _colors["gray53"] = RGBColor(135, 135, 135);
  _colors["gray54"] = RGBColor(138, 138, 138);
  _colors["gray55"] = RGBColor(140, 140, 140);
  _colors["gray56"] = RGBColor(143, 143, 143);
  _colors["gray57"] = RGBColor(145, 145, 145);
  _colors["gray58"] = RGBColor(148, 148, 148);
  _colors["gray59"] = RGBColor(150, 150, 150);
  _colors["gray60"] = RGBColor(153, 153, 153);
  _colors["gray61"] = RGBColor(156, 156, 156);
  _colors["gray62"] = RGBColor(158, 158, 158);
  _colors["gray63"] = RGBColor(161, 161, 161);
  _colors["gray64"] = RGBColor(163, 163, 163);
  _colors["gray65"] = RGBColor(166, 166, 166);
  _colors["gray66"] = RGBColor(168, 168, 168);
  _colors["gray67"] = RGBColor(171, 171, 171);
  _colors["gray68"] = RGBColor(173, 173, 173);
  _colors["gray69"] = RGBColor(176, 176, 176);
  _colors["gray70"] = RGBColor(179, 179, 179);
  _colors["gray71"] = RGBColor(181, 181, 181);
  _colors["gray72"] = RGBColor(184, 184, 184);
  _colors["gray73"] = RGBColor(186, 186, 186);
  _colors["gray74"] = RGBColor(189, 189, 189);
  _colors["gray75"] = RGBColor(191, 191, 191);
  _colors["gray76"] = RGBColor(194, 194, 194);
  _colors["gray77"] = RGBColor(196, 196, 196);
  _colors["gray78"] = RGBColor(199, 199, 199);
  _colors["gray79"] = RGBColor(201, 201, 201);
  _colors["gray80"] = RGBColor(204, 204, 204);
  _colors["gray81"] = RGBColor(207, 207, 207);
  _colors["gray82"] = RGBColor(209, 209, 209);
  _colors["gray83"] = RGBColor(212, 212, 212);
  _colors["gray84"] = RGBColor(214, 214, 214);
  _colors["gray85"] = RGBColor(217, 217, 217);
  _colors["gray86"] = RGBColor(219, 219, 219);
  _colors["gray87"] = RGBColor(222, 222, 222);
  _colors["gray88"] = RGBColor(224, 224, 224);
  _colors["gray89"] = RGBColor(227, 227, 227);
  _colors["gray90"] = RGBColor(229, 229, 229);
  _colors["gray91"] = RGBColor(232, 232, 232);
  _colors["gray92"] = RGBColor(235, 235, 235);
  _colors["gray93"] = RGBColor(237, 237, 237);
  _colors["gray94"] = RGBColor(240, 240, 240);
  _colors["gray95"] = RGBColor(242, 242, 242);
  _colors["gray96"] = RGBColor(245, 245, 245);
  _colors["gray97"] = RGBColor(247, 247, 247);
  _colors["gray98"] = RGBColor(250, 250, 250);
  _colors["gray99"] = RGBColor(252, 252, 252);
  _colors["gray100"] = RGBColor(255, 255, 255);
  _colors["green"] = RGBColor(0, 255, 0);
  _colors["green1"] = RGBColor(0, 255, 0);
  _colors["green2"] = RGBColor(0, 238, 0);
  _colors["green3"] = RGBColor(0, 205, 0);
  _colors["green4"] = RGBColor(0, 139, 0);
  _colors["greenyellow"] = RGBColor(173, 255, 47);
  _colors["grey"] = RGBColor(190, 190, 190);
  _colors["grey0"] = RGBColor(0, 0, 0);
  _colors["grey1"] = RGBColor(3, 3, 3);
  _colors["grey2"] = RGBColor(5, 5, 5);
  _colors["grey3"] = RGBColor(8, 8, 8);
  _colors["grey4"] = RGBColor(10, 10, 10);
  _colors["grey5"] = RGBColor(13, 13, 13);
  _colors["grey6"] = RGBColor(15, 15, 15);
  _colors["grey7"] = RGBColor(18, 18, 18);
  _colors["grey8"] = RGBColor(20, 20, 20);
  _colors["grey9"] = RGBColor(23, 23, 23);
  _colors["grey10"] = RGBColor(26, 26, 26);
  _colors["grey11"] = RGBColor(28, 28, 28);
  _colors["grey12"] = RGBColor(31, 31, 31);
  _colors["grey13"] = RGBColor(33, 33, 33);
  _colors["grey14"] = RGBColor(36, 36, 36);
  _colors["grey15"] = RGBColor(38, 38, 38);
  _colors["grey16"] = RGBColor(41, 41, 41);
  _colors["grey17"] = RGBColor(43, 43, 43);
  _colors["grey18"] = RGBColor(46, 46, 46);
  _colors["grey19"] = RGBColor(48, 48, 48);
  _colors["grey20"] = RGBColor(51, 51, 51);
  _colors["grey21"] = RGBColor(54, 54, 54);
  _colors["grey22"] = RGBColor(56, 56, 56);
  _colors["grey23"] = RGBColor(59, 59, 59);
  _colors["grey24"] = RGBColor(61, 61, 61);
  _colors["grey25"] = RGBColor(64, 64, 64);
  _colors["grey26"] = RGBColor(66, 66, 66);
  _colors["grey27"] = RGBColor(69, 69, 69);
  _colors["grey28"] = RGBColor(71, 71, 71);
  _colors["grey29"] = RGBColor(74, 74, 74);
  _colors["grey30"] = RGBColor(77, 77, 77);
  _colors["grey31"] = RGBColor(79, 79, 79);
  _colors["grey32"] = RGBColor(82, 82, 82);
  _colors["grey33"] = RGBColor(84, 84, 84);
  _colors["grey34"] = RGBColor(87, 87, 87);
  _colors["grey35"] = RGBColor(89, 89, 89);
  _colors["grey36"] = RGBColor(92, 92, 92);
  _colors["grey37"] = RGBColor(94, 94, 94);
  _colors["grey38"] = RGBColor(97, 97, 97);
  _colors["grey39"] = RGBColor(99, 99, 99);
  _colors["grey40"] = RGBColor(102, 102, 102);
  _colors["grey41"] = RGBColor(105, 105, 105);
  _colors["grey42"] = RGBColor(107, 107, 107);
  _colors["grey43"] = RGBColor(110, 110, 110);
  _colors["grey44"] = RGBColor(112, 112, 112);
  _colors["grey45"] = RGBColor(115, 115, 115);
  _colors["grey46"] = RGBColor(117, 117, 117);
  _colors["grey47"] = RGBColor(120, 120, 120);
  _colors["grey48"] = RGBColor(122, 122, 122);
  _colors["grey49"] = RGBColor(125, 125, 125);
  _colors["grey50"] = RGBColor(127, 127, 127);
  _colors["grey51"] = RGBColor(130, 130, 130);
  _colors["grey52"] = RGBColor(133, 133, 133);
  _colors["grey53"] = RGBColor(135, 135, 135);
  _colors["grey54"] = RGBColor(138, 138, 138);
  _colors["grey55"] = RGBColor(140, 140, 140);
  _colors["grey56"] = RGBColor(143, 143, 143);
  _colors["grey57"] = RGBColor(145, 145, 145);
  _colors["grey58"] = RGBColor(148, 148, 148);
  _colors["grey59"] = RGBColor(150, 150, 150);
  _colors["grey60"] = RGBColor(153, 153, 153);
  _colors["grey61"] = RGBColor(156, 156, 156);
  _colors["grey62"] = RGBColor(158, 158, 158);
  _colors["grey63"] = RGBColor(161, 161, 161);
  _colors["grey64"] = RGBColor(163, 163, 163);
  _colors["grey65"] = RGBColor(166, 166, 166);
  _colors["grey66"] = RGBColor(168, 168, 168);
  _colors["grey67"] = RGBColor(171, 171, 171);
  _colors["grey68"] = RGBColor(173, 173, 173);
  _colors["grey69"] = RGBColor(176, 176, 176);
  _colors["grey70"] = RGBColor(179, 179, 179);
  _colors["grey71"] = RGBColor(181, 181, 181);
  _colors["grey72"] = RGBColor(184, 184, 184);
  _colors["grey73"] = RGBColor(186, 186, 186);
  _colors["grey74"] = RGBColor(189, 189, 189);
  _colors["grey75"] = RGBColor(191, 191, 191);
  _colors["grey76"] = RGBColor(194, 194, 194);
  _colors["grey77"] = RGBColor(196, 196, 196);
  _colors["grey78"] = RGBColor(199, 199, 199);
  _colors["grey79"] = RGBColor(201, 201, 201);
  _colors["grey80"] = RGBColor(204, 204, 204);
  _colors["grey81"] = RGBColor(207, 207, 207);
  _colors["grey82"] = RGBColor(209, 209, 209);
  _colors["grey83"] = RGBColor(212, 212, 212);
  _colors["grey84"] = RGBColor(214, 214, 214);
  _colors["grey85"] = RGBColor(217, 217, 217);
  _colors["grey86"] = RGBColor(219, 219, 219);
  _colors["grey87"] = RGBColor(222, 222, 222);
  _colors["grey88"] = RGBColor(224, 224, 224);
  _colors["grey89"] = RGBColor(227, 227, 227);
  _colors["grey90"] = RGBColor(229, 229, 229);
  _colors["grey91"] = RGBColor(232, 232, 232);
  _colors["grey92"] = RGBColor(235, 235, 235);
  _colors["grey93"] = RGBColor(237, 237, 237);
  _colors["grey94"] = RGBColor(240, 240, 240);
  _colors["grey95"] = RGBColor(242, 242, 242);
  _colors["grey96"] = RGBColor(245, 245, 245);
  _colors["grey97"] = RGBColor(247, 247, 247);
  _colors["grey98"] = RGBColor(250, 250, 250);
  _colors["grey99"] = RGBColor(252, 252, 252);
  _colors["grey100"] = RGBColor(255, 255, 255);
  _colors["honeydew"] = RGBColor(240, 255, 240);
  _colors["honeydew1"] = RGBColor(240, 255, 240);
  _colors["honeydew2"] = RGBColor(224, 238, 224);
  _colors["honeydew3"] = RGBColor(193, 205, 193);
  _colors["honeydew4"] = RGBColor(131, 139, 131);
  _colors["hotpink"] = RGBColor(255, 105, 180);
  _colors["hotpink1"] = RGBColor(255, 110, 180);
  _colors["hotpink2"] = RGBColor(238, 106, 167);
  _colors["hotpink3"] = RGBColor(205, 96, 144);
  _colors["hotpink4"] = RGBColor(139, 58, 98);
  _colors["indianred"] = RGBColor(205, 92, 92);
  _colors["indianred1"] = RGBColor(255, 106, 106);
  _colors["indianred2"] = RGBColor(238, 99, 99);
  _colors["indianred3"] = RGBColor(205, 85, 85);
  _colors["indianred4"] = RGBColor(139, 58, 58);
  _colors["ivory"] = RGBColor(255, 255, 240);
  _colors["ivory1"] = RGBColor(255, 255, 240);
  _colors["ivory2"] = RGBColor(238, 238, 224);
  _colors["ivory3"] = RGBColor(205, 205, 193);
  _colors["ivory4"] = RGBColor(139, 139, 131);
  _colors["khaki"] = RGBColor(240, 230, 140);
  _colors["khaki1"] = RGBColor(255, 246, 143);
  _colors["khaki2"] = RGBColor(238, 230, 133);
  _colors["khaki3"] = RGBColor(205, 198, 115);
  _colors["khaki4"] = RGBColor(139, 134, 78);
  _colors["lavender"] = RGBColor(230, 230, 250);
  _colors["lavenderblush"] = RGBColor(255, 240, 245);
  _colors["lavenderblush1"] = RGBColor(255, 240, 245);
  _colors["lavenderblush2"] = RGBColor(238, 224, 229);
  _colors["lavenderblush3"] = RGBColor(205, 193, 197);
  _colors["lavenderblush4"] = RGBColor(139, 131, 134);
  _colors["lawngreen"] = RGBColor(124, 252, 0);
  _colors["lemonchiffon"] = RGBColor(255, 250, 205);
  _colors["lemonchiffon1"] = RGBColor(255, 250, 205);
  _colors["lemonchiffon2"] = RGBColor(238, 233, 191);
  _colors["lemonchiffon3"] = RGBColor(205, 201, 165);
  _colors["lemonchiffon4"] = RGBColor(139, 137, 112);
  _colors["lightblue"] = RGBColor(173, 216, 230);
  _colors["lightblue1"] = RGBColor(191, 239, 255);
  _colors["lightblue2"] = RGBColor(178, 223, 238);
  _colors["lightblue3"] = RGBColor(154, 192, 205);
  _colors["lightblue4"] = RGBColor(104, 131, 139);
  _colors["lightcoral"] = RGBColor(240, 128, 128);
  _colors["lightcyan"] = RGBColor(224, 255, 255);
  _colors["lightcyan1"] = RGBColor(224, 255, 255);
  _colors["lightcyan2"] = RGBColor(209, 238, 238);
  _colors["lightcyan3"] = RGBColor(180, 205, 205);
  _colors["lightcyan4"] = RGBColor(122, 139, 139);
  _colors["lightgoldenrod"] = RGBColor(238, 221, 130);
  _colors["lightgoldenrod1"] = RGBColor(255, 236, 139);
  _colors["lightgoldenrod2"] = RGBColor(238, 220, 130);
  _colors["lightgoldenrod3"] = RGBColor(205, 190, 112);
  _colors["lightgoldenrod4"] = RGBColor(139, 129, 76);
  _colors["lightgoldenrodyellow"] = RGBColor(250, 250, 210);
  _colors["lightgray"] = RGBColor(211, 211, 211);
  _colors["lightgreen"] = RGBColor(144, 238, 144);
  _colors["lightgrey"] = RGBColor(211, 211, 211);
  _colors["lightpink"] = RGBColor(255, 182, 193);
  _colors["lightpink1"] = RGBColor(255, 174, 185);
  _colors["lightpink2"] = RGBColor(238, 162, 173);
  _colors["lightpink3"] = RGBColor(205, 140, 149);
  _colors["lightpink4"] = RGBColor(139, 95, 101);
  _colors["lightsalmon"] = RGBColor(255, 160, 122);
  _colors["lightsalmon1"] = RGBColor(255, 160, 122);
  _colors["lightsalmon2"] = RGBColor(238, 149, 114);
  _colors["lightsalmon3"] = RGBColor(205, 129, 98);
  _colors["lightsalmon4"] = RGBColor(139, 87, 66);
  _colors["lightseagreen"] = RGBColor(32, 178, 170);
  _colors["lightskyblue"] = RGBColor(135, 206, 250);
  _colors["lightskyblue1"] = RGBColor(176, 226, 255);
  _colors["lightskyblue2"] = RGBColor(164, 211, 238);
  _colors["lightskyblue3"] = RGBColor(141, 182, 205);
  _colors["lightskyblue4"] = RGBColor(96, 123, 139);
  _colors["lightslateblue"] = RGBColor(132, 112, 255);
  _colors["lightslategray"] = RGBColor(119, 136, 153);
  _colors["lightslategrey"] = RGBColor(119, 136, 153);
  _colors["lightsteelblue"] = RGBColor(176, 196, 222);
  _colors["lightsteelblue1"] = RGBColor(202, 225, 255);
  _colors["lightsteelblue2"] = RGBColor(188, 210, 238);
  _colors["lightsteelblue3"] = RGBColor(162, 181, 205);
  _colors["lightsteelblue4"] = RGBColor(110, 123, 139);
  _colors["lightyellow"] = RGBColor(255, 255, 224);
  _colors["lightyellow1"] = RGBColor(255, 255, 224);
  _colors["lightyellow2"] = RGBColor(238, 238, 209);
  _colors["lightyellow3"] = RGBColor(205, 205, 180);
  _colors["lightyellow4"] = RGBColor(139, 139, 122);
  _colors["limegreen"] = RGBColor(50, 205, 50);
  _colors["linen"] = RGBColor(250, 240, 230);
  _colors["magenta"] = RGBColor(255, 0, 255);
  _colors["magenta1"] = RGBColor(255, 0, 255);
  _colors["magenta2"] = RGBColor(238, 0, 238);
  _colors["magenta3"] = RGBColor(205, 0, 205);
  _colors["magenta4"] = RGBColor(139, 0, 139);
  _colors["maroon"] = RGBColor(176, 48, 96);
  _colors["maroon1"] = RGBColor(255, 52, 179);
  _colors["maroon2"] = RGBColor(238, 48, 167);
  _colors["maroon3"] = RGBColor(205, 41, 144);
  _colors["maroon4"] = RGBColor(139, 28, 98);
  _colors["mediumaquamarine"] = RGBColor(102, 205, 170);
  _colors["mediumblue"] = RGBColor(0, 0, 205);
  _colors["mediumorchid"] = RGBColor(186, 85, 211);
  _colors["mediumorchid1"] = RGBColor(224, 102, 255);
  _colors["mediumorchid2"] = RGBColor(209, 95, 238);
  _colors["mediumorchid3"] = RGBColor(180, 82, 205);
  _colors["mediumorchid4"] = RGBColor(122, 55, 139);
  _colors["mediumpurple"] = RGBColor(147, 112, 219);
  _colors["mediumpurple1"] = RGBColor(171, 130, 255);
  _colors["mediumpurple2"] = RGBColor(159, 121, 238);
  _colors["mediumpurple3"] = RGBColor(137, 104, 205);
  _colors["mediumpurple4"] = RGBColor(93, 71, 139);
  _colors["mediumseagreen"] = RGBColor(60, 179, 113);
  _colors["mediumslateblue"] = RGBColor(123, 104, 238);
  _colors["mediumspringgreen"] = RGBColor(0, 250, 154);
  _colors["mediumturquoise"] = RGBColor(72, 209, 204);
  _colors["mediumvioletred"] = RGBColor(199, 21, 133);
  _colors["midnightblue"] = RGBColor(25, 25, 112);
  _colors["mintcream"] = RGBColor(245, 255, 250);
  _colors["mistyrose"] = RGBColor(255, 228, 225);
  _colors["mistyrose1"] = RGBColor(255, 228, 225);
  _colors["mistyrose2"] = RGBColor(238, 213, 210);
  _colors["mistyrose3"] = RGBColor(205, 183, 181);
  _colors["mistyrose4"] = RGBColor(139, 125, 123);
  _colors["moccasin"] = RGBColor(255, 228, 181);
  _colors["navajowhite"] = RGBColor(255, 222, 173);
  _colors["navajowhite1"] = RGBColor(255, 222, 173);
  _colors["navajowhite2"] = RGBColor(238, 207, 161);
  _colors["navajowhite3"] = RGBColor(205, 179, 139);
  _colors["navajowhite4"] = RGBColor(139, 121, 94);
  _colors["navy"] = RGBColor(0, 0, 128);
  _colors["navyblue"] = RGBColor(0, 0, 128);
  _colors["oldlace"] = RGBColor(253, 245, 230);
  _colors["olivedrab"] = RGBColor(107, 142, 35);
  _colors["olivedrab1"] = RGBColor(192, 255, 62);
  _colors["olivedrab2"] = RGBColor(179, 238, 58);
  _colors["olivedrab3"] = RGBColor(154, 205, 50);
  _colors["olivedrab4"] = RGBColor(105, 139, 34);
  _colors["orange"] = RGBColor(255, 165, 0);
  _colors["orange1"] = RGBColor(255, 165, 0);
  _colors["orange2"] = RGBColor(238, 154, 0);
  _colors["orange3"] = RGBColor(205, 133, 0);
  _colors["orange4"] = RGBColor(139, 90, 0);
  _colors["orangered"] = RGBColor(255, 69, 0);
  _colors["orangered1"] = RGBColor(255, 69, 0);
  _colors["orangered2"] = RGBColor(238, 64, 0);
  _colors["orangered3"] = RGBColor(205, 55, 0);
  _colors["orangered4"] = RGBColor(139, 37, 0);
  _colors["orchid"] = RGBColor(218, 112, 214);
  _colors["orchid1"] = RGBColor(255, 131, 250);
  _colors["orchid2"] = RGBColor(238, 122, 233);
  _colors["orchid3"] = RGBColor(205, 105, 201);
  _colors["orchid4"] = RGBColor(139, 71, 137);
  _colors["palegoldenrod"] = RGBColor(238, 232, 170);
  _colors["palegreen"] = RGBColor(152, 251, 152);
  _colors["palegreen1"] = RGBColor(154, 255, 154);
  _colors["palegreen2"] = RGBColor(144, 238, 144);
  _colors["palegreen3"] = RGBColor(124, 205, 124);
  _colors["palegreen4"] = RGBColor(84, 139, 84);
  _colors["paleturquoise"] = RGBColor(175, 238, 238);
  _colors["paleturquoise1"] = RGBColor(187, 255, 255);
  _colors["paleturquoise2"] = RGBColor(174, 238, 238);
  _colors["paleturquoise3"] = RGBColor(150, 205, 205);
  _colors["paleturquoise4"] = RGBColor(102, 139, 139);
  _colors["palevioletred"] = RGBColor(219, 112, 147);
  _colors["palevioletred1"] = RGBColor(255, 130, 171);
  _colors["palevioletred2"] = RGBColor(238, 121, 159);
  _colors["palevioletred3"] = RGBColor(205, 104, 137);
  _colors["palevioletred4"] = RGBColor(139, 71, 93);
  _colors["papayawhip"] = RGBColor(255, 239, 213);
  _colors["peachpuff"] = RGBColor(255, 218, 185);
  _colors["peachpuff1"] = RGBColor(255, 218, 185);
  _colors["peachpuff2"] = RGBColor(238, 203, 173);
  _colors["peachpuff3"] = RGBColor(205, 175, 149);
  _colors["peachpuff4"] = RGBColor(139, 119, 101);
  _colors["peru"] = RGBColor(205, 133, 63);
  _colors["pink"] = RGBColor(255, 192, 203);
  _colors["pink1"] = RGBColor(255, 181, 197);
  _colors["pink2"] = RGBColor(238, 169, 184);
  _colors["pink3"] = RGBColor(205, 145, 158);
  _colors["pink4"] = RGBColor(139, 99, 108);
  _colors["plum"] = RGBColor(221, 160, 221);
  _colors["plum1"] = RGBColor(255, 187, 255);
  _colors["plum2"] = RGBColor(238, 174, 238);
  _colors["plum3"] = RGBColor(205, 150, 205);
  _colors["plum4"] = RGBColor(139, 102, 139);
  _colors["powderblue"] = RGBColor(176, 224, 230);
  _colors["purple"] = RGBColor(160, 32, 240);
  _colors["purple1"] = RGBColor(155, 48, 255);
  _colors["purple2"] = RGBColor(145, 44, 238);
  _colors["purple3"] = RGBColor(125, 38, 205);
  _colors["purple4"] = RGBColor(85, 26, 139);
  _colors["red"] = RGBColor(255, 0, 0);
  _colors["red1"] = RGBColor(255, 0, 0);
  _colors["red2"] = RGBColor(238, 0, 0);
  _colors["red3"] = RGBColor(205, 0, 0);
  _colors["red4"] = RGBColor(139, 0, 0);
  _colors["rosybrown"] = RGBColor(188, 143, 143);
  _colors["rosybrown1"] = RGBColor(255, 193, 193);
  _colors["rosybrown2"] = RGBColor(238, 180, 180);
  _colors["rosybrown3"] = RGBColor(205, 155, 155);
  _colors["rosybrown4"] = RGBColor(139, 105, 105);
  _colors["royalblue"] = RGBColor(65, 105, 225);
  _colors["royalblue1"] = RGBColor(72, 118, 255);
  _colors["royalblue2"] = RGBColor(67, 110, 238);
  _colors["royalblue3"] = RGBColor(58, 95, 205);
  _colors["royalblue4"] = RGBColor(39, 64, 139);
  _colors["saddlebrown"] = RGBColor(139, 69, 19);
  _colors["salmon"] = RGBColor(250, 128, 114);
  _colors["salmon1"] = RGBColor(255, 140, 105);
  _colors["salmon2"] = RGBColor(238, 130, 98);
  _colors["salmon3"] = RGBColor(205, 112, 84);
  _colors["salmon4"] = RGBColor(139, 76, 57);
  _colors["sandybrown"] = RGBColor(244, 164, 96);
  _colors["seagreen"] = RGBColor(46, 139, 87);
  _colors["seagreen1"] = RGBColor(84, 255, 159);
  _colors["seagreen2"] = RGBColor(78, 238, 148);
  _colors["seagreen3"] = RGBColor(67, 205, 128);
  _colors["seagreen4"] = RGBColor(46, 139, 87);
  _colors["seashell"] = RGBColor(255, 245, 238);
  _colors["seashell1"] = RGBColor(255, 245, 238);
  _colors["seashell2"] = RGBColor(238, 229, 222);
  _colors["seashell3"] = RGBColor(205, 197, 191);
  _colors["seashell4"] = RGBColor(139, 134, 130);
  _colors["sienna"] = RGBColor(160, 82, 45);
  _colors["sienna1"] = RGBColor(255, 130, 71);
  _colors["sienna2"] = RGBColor(238, 121, 66);
  _colors["sienna3"] = RGBColor(205, 104, 57);
  _colors["sienna4"] = RGBColor(139, 71, 38);
  _colors["skyblue"] = RGBColor(135, 206, 235);
  _colors["skyblue1"] = RGBColor(135, 206, 255);
  _colors["skyblue2"] = RGBColor(126, 192, 238);
  _colors["skyblue3"] = RGBColor(108, 166, 205);
  _colors["skyblue4"] = RGBColor(74, 112, 139);
  _colors["slateblue"] = RGBColor(106, 90, 205);
  _colors["slateblue1"] = RGBColor(131, 111, 255);
  _colors["slateblue2"] = RGBColor(122, 103, 238);
  _colors["slateblue3"] = RGBColor(105, 89, 205);
  _colors["slateblue4"] = RGBColor(71, 60, 139);
  _colors["slategray"] = RGBColor(112, 128, 144);
  _colors["slategray1"] = RGBColor(198, 226, 255);
  _colors["slategray2"] = RGBColor(185, 211, 238);
  _colors["slategray3"] = RGBColor(159, 182, 205);
  _colors["slategray4"] = RGBColor(108, 123, 139);
  _colors["slategrey"] = RGBColor(112, 128, 144);
  _colors["snow"] = RGBColor(255, 250, 250);
  _colors["snow1"] = RGBColor(255, 250, 250);
  _colors["snow2"] = RGBColor(238, 233, 233);
  _colors["snow3"] = RGBColor(205, 201, 201);
  _colors["snow4"] = RGBColor(139, 137, 137);
  _colors["springgreen"] = RGBColor(0, 255, 127);
  _colors["springgreen1"] = RGBColor(0, 255, 127);
  _colors["springgreen2"] = RGBColor(0, 238, 118);
  _colors["springgreen3"] = RGBColor(0, 205, 102);
  _colors["springgreen4"] = RGBColor(0, 139, 69);
  _colors["steelblue"] = RGBColor(70, 130, 180);
  _colors["steelblue1"] = RGBColor(99, 184, 255);
  _colors["steelblue2"] = RGBColor(92, 172, 238);
  _colors["steelblue3"] = RGBColor(79, 148, 205);
  _colors["steelblue4"] = RGBColor(54, 100, 139);
  _colors["tan"] = RGBColor(210, 180, 140);
  _colors["tan1"] = RGBColor(255, 165, 79);
  _colors["tan2"] = RGBColor(238, 154, 73);
  _colors["tan3"] = RGBColor(205, 133, 63);
  _colors["tan4"] = RGBColor(139, 90, 43);
  _colors["thistle"] = RGBColor(216, 191, 216);
  _colors["thistle1"] = RGBColor(255, 225, 255);
  _colors["thistle2"] = RGBColor(238, 210, 238);
  _colors["thistle3"] = RGBColor(205, 181, 205);
  _colors["thistle4"] = RGBColor(139, 123, 139);
  _colors["tomato"] = RGBColor(255, 99, 71);
  _colors["tomato1"] = RGBColor(255, 99, 71);
  _colors["tomato2"] = RGBColor(238, 92, 66);
  _colors["tomato3"] = RGBColor(205, 79, 57);
  _colors["tomato4"] = RGBColor(139, 54, 38);
  _colors["turquoise"] = RGBColor(64, 224, 208);
  _colors["turquoise1"] = RGBColor(0, 245, 255);
  _colors["turquoise2"] = RGBColor(0, 229, 238);
  _colors["turquoise3"] = RGBColor(0, 197, 205);
  _colors["turquoise4"] = RGBColor(0, 134, 139);
  _colors["violet"] = RGBColor(238, 130, 238);
  _colors["violetred"] = RGBColor(208, 32, 144);
  _colors["violetred1"] = RGBColor(255, 62, 150);
  _colors["violetred2"] = RGBColor(238, 58, 140);
  _colors["violetred3"] = RGBColor(205, 50, 120);
  _colors["violetred4"] = RGBColor(139, 34, 82);
  _colors["wheat"] = RGBColor(245, 222, 179);
  _colors["wheat1"] = RGBColor(255, 231, 186);
  _colors["wheat2"] = RGBColor(238, 216, 174);
  _colors["wheat3"] = RGBColor(205, 186, 150);
  _colors["wheat4"] = RGBColor(139, 126, 102);
  _colors["whitesmoke"] = RGBColor(245, 245, 245);
  _colors["yellow"] = RGBColor(255, 255, 0);
  _colors["yellow1"] = RGBColor(255, 255, 0);
  _colors["yellow2"] = RGBColor(238, 238, 0);
  _colors["yellow3"] = RGBColor(205, 205, 0);
  _colors["yellow4"] = RGBColor(139, 139, 0);
  _colors["yellowgreen"] = RGBColor(154, 205, 50);
}

