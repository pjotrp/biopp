//
// File: Newick.h
// Created by: Julien Dutheil
// Created on: Thu Oct 23 15:35:03 2003
//

/*
Copyright or � or Copr. CNRS, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

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

#ifndef _NEWICK_H_
#define _NEWICK_H_

#include "IOTree.h"
#include "TreeTemplate.h"

namespace bpp
{

/**
 * @brief The so-called 'newick' parenthetic format.
 *
 * Branch lengths and bootstraps are supported:
 *
 * ex:
 * <code>
 * ((A:0.1, B:0.15)90:0.2, C:0.27);
 * </code>
 *
 * Code example:
 * @code
 * #include <Phyl/Newick.h>
 * #include <Phyl/Tree.h>
 * 
 * Newick * newickReader = new Newick(false); //No comment allowed!
 * try {
 * 	Tree * tree = newickReader->read("MyTestTree.dnd"); // Tree in file MyTestTree.dnd
 * 	cout << "Tree has " << tree->getNumberOfLeaves() << " leaves." << endl;
 * } catch (Exception e) {
 *	cout << "Error when reading tree." << endl;
 * }
 * delete tree;
 * delete newickReader;
 * @endcode
 *
 * Bootstrap values are stored as node properties, as Number<double> objects and with the tag TreeTools::BOOTSTRAP.
 *
 * This is also possible to read a "tagged" tree, where additional info is provided in place of bootstrap values:
 * ((A,B)N2,(C,D)N3)N1;
 * This is achieved by calling the enableExtendedBootstrapProperty method, and porviding a property name to use.
 * The additional information will be stored at each node as a property, in a String object.
 * The disableExtendedBootstrapProperty method restores the default behavior.
 */
class Newick:
  public virtual AbstractITree,
  public virtual AbstractOTree,
  public virtual AbstractIMultiTree,
  public virtual AbstractOMultiTree
{
	protected:
		bool _allowComments;
    bool _writeId;
    bool _useBootstrap;
    string _bootstrapPropertyName;
	
	public:
		
		/**
		 * @brief Build a new Newick reader/writer.
		 *
		 * Some newick format allow comments between hooks ('[' ']').
		 * 
		 * @param allowComments Tell if comments between [] are allowed in file.
		 * @param writeId       If true, nodes ids will be written in place of bootstrap values.
		 */
		Newick(bool allowComments = false, bool writeId = false):
      _allowComments(allowComments),
      _writeId(writeId),
      _useBootstrap(true),
      _bootstrapPropertyName(TreeTools::BOOTSTRAP) {}

		virtual ~Newick() {}
	
	public:

    void enableExtendedBootstrapProperty(const string & propertyName)
    {
      _useBootstrap = false;
      _bootstrapPropertyName = propertyName;
    }
    void disableExtendedBootstrapProperty()
    {
      _useBootstrap = true;
      _bootstrapPropertyName = TreeTools::BOOTSTRAP;
    }

		/**
		 * @name The IOTree interface
		 *
		 * @{
		 */
		const string getFormatName() const;
		const string getFormatDescription() const;
		/* @} */

		/**
		 * @name The ITree interface
		 *
		 * @{
		 */
		
#if defined(NO_VIRTUAL_COV)
		Tree * read(const string & path) const throw (Exception)
		{
			return AbstractITree::read(path);
		}
#else
		TreeTemplate<Node> * read(const string & path) const throw (Exception)
		{
			return dynamic_cast<TreeTemplate<Node> *>(AbstractITree::read(path));
		}
#endif
		
#if defined(NO_VIRTUAL_COV)
		Tree *
#else
		TreeTemplate<Node> * 
#endif
		read(istream & in) const throw (Exception);
		/** @} */

		/**
		 * @name The OTree interface
		 *
		 * @{
		 */
		void write(const Tree & tree, const string & path, bool overwrite = true) const throw (Exception)
		{
			AbstractOTree::write(tree, path, overwrite);
		}
		void write(const Tree & tree, ostream & out) const throw (Exception);
		/** @} */

		/**
		 * @name The IMultiTree interface
		 *
		 * @{
		 */
		void read(const string & path, vector<Tree *> & trees) const throw (Exception)
		{
			AbstractIMultiTree::read(path, trees);
		}
		void read(istream & in, vector<Tree *> & trees) const throw (Exception);
    /**@}*/

		/**
		 * @name The OMultiTree interface
		 *
		 * @{
		 */
		void write(const vector<Tree *> & trees, const string & path, bool overwrite = true) const throw (Exception)
		{
			AbstractOMultiTree::write(trees, path, overwrite);
		}
		void write(const vector<Tree *> & trees, ostream & out) const throw (Exception);
		/** @} */
};

} //end of namespace bpp.

#endif	//_NEWICK_H_

