//
// File: DataTable.h
// Created by: Julien Dutheil
// Created on: Aug 2005
//

/*
Copyright or � or Copr. CNRS, (November 17, 2004)

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

#ifndef _DataTable_H_
#define _DataTable_H_

#include "VectorTools.h"
#include "DataTableExceptions.h"

// From Utils:
#include <Utils/Exceptions.h>
#include <Utils/TextTools.h>
#include <Utils/Clonable.h>

// From the STL:
#include <string>
#include <map>
using namespace std;

namespace bpp
{

/**
 * @brief This class corresponds to a 'dataset', <i>i.e.</i> a table with data by rows
 * and variable by columns.
 *
 * Data are stored as string objects, by column.
 * A DataTable object is hence similar to a ColMatrix<string>.object.
 * (NB: actually, ColMatrix does not exist yet...)
 */
class DataTable:
  public Clonable
{
  
  protected:
    unsigned int nRow_, nCol_;
    vector< vector<string> > data_;
    vector<string> * rowNames_;
    vector<string> * colNames_;

  public:
    
    /**
     * @brief Build a new void DataTable object with nRow rows and nCol columns.
     *
     * @param nRow The number of rows of the DataTable.
     * @param nCol The number of columns of the DataTable.
     */
    DataTable(unsigned int nRow, unsigned int nCol);
    
    /**
     * @brief Build a new void DataTable object with nCol columns.
     *
     * @param nCol The number of columns of the DataTable.
     */
    DataTable(unsigned int nCol);

    /**
     * @brief Build a new void DataTable object with named columns.
     *
     * @param colNames The names of the columns of the DataTable.
     * @throw DuplicatedTableColumnNameException If colnames contains identical names.
     */
    DataTable(const vector<string> & colNames) throw (DuplicatedTableColumnNameException);

    DataTable(const DataTable & table);

    DataTable & operator=(const DataTable & table);

    DataTable * clone() const { return new DataTable(*this); }

    virtual ~DataTable();

  public:

    /**
     * @return The element at a given position.
     * @param rowIndex Row number.
     * @param colIndex Column number.
     * @throw IndexOutOfBoundsException If one of the index is greater or equal to the corresponding number of columns/rows. 
     */
    string & operator()(unsigned int rowIndex, unsigned int colIndex) throw (IndexOutOfBoundsException);
    
    /**
     * @return The element at a given position.
     * @param rowIndex Row number.
     * @param colIndex Column number.
     * @throw IndexOutOfBoundsException If one of the index is greater or equal to the corresponding number of columns/rows. 
     */
    const string & operator()(unsigned int rowIndex, unsigned int colIndex) const throw (IndexOutOfBoundsException);
    
    /**
     * @return The element at a given position.
     * @param rowName Row name.
     * @param colName Column name.
     * @throw NoTableRowNamesException If the table does not have names associated to rows. 
     * @throw NoTableColumnNamesException If the table does not have names associated to columns. 
     * @throw TableNameNotFoundException If one of rowName or colName do not match existing names. 
     */
    string & operator()(const string & rowName, const string & colName)
            throw (NoTableRowNamesException, NoTableColumnNamesException, TableNameNotFoundException);
    
    /**
     * @return The element at a given position.
     * @param rowName Row name.
     * @param colName Column name.
     * @throw NoTableRowNamesException If the table does not have names associated to rows. 
     * @throw NoTableColumnNamesException If the table does not have names associated to columns. 
     * @throw TableNameNotFoundException If one of rowName or colName do not match existing names. 
     */
    const string & operator()(const string & rowName, const string & colName) const
            throw (NoTableRowNamesException, NoTableColumnNamesException, TableNameNotFoundException);
    
    /**
     * @return The element at a given position.
     * @param rowName Row name.
     * @param colIndex Column number.
     * @throw NoTableRowNamesException If the table does not have names associated to rows. 
     * @throw IndexOutOfBoundsException If the index is greater or equal to the number of columns. 
     * @throw TableNameNotFoundException If rowName do not match existing names. 
     */
    string & operator()(const string & rowName, unsigned int colIndex)
            throw (NoTableRowNamesException, TableNameNotFoundException, IndexOutOfBoundsException);
    
    /**
     * @return The element at a given position.
     * @param rowName Row name.
     * @param colIndex Column number.
     * @throw NoTableRowNamesException If the table does not have names associated to rows. 
     * @throw IndexOutOfBoundsException If the index is greater or equal to the number of columns. 
     * @throw TableNameNotFoundException If rowName do not match existing names. 
     */
    const string & operator()(const string & rowName, unsigned int colIndex) const
            throw (NoTableRowNamesException, TableNameNotFoundException, IndexOutOfBoundsException);

    /**
     * @return The element at a given position.
     * @param rowIndex Row number.
     * @param colName Column name.
     * @throw IndexOutOfBoundsException If the index is greater or equal to the number of rows. 
     * @throw NoTableColumnNamesException If the table does not have names associated to columns. 
     * @throw TableNameNotFoundException If colName do not match existing names. 
     */
    string & operator()(unsigned int rowIndex, const string & colName)
            throw (IndexOutOfBoundsException, NoTableColumnNamesException, TableNameNotFoundException);
    
    /**
     * @return The element at a given position.
     * @param rowIndex Row number.
     * @param colName Column name.
     * @throw IndexOutOfBoundsException If the index is greater or equal to the number of rows. 
     * @throw NoTableColumnNamesException If the table does not have names associated to columns. 
     * @throw TableNameNotFoundException If colName do not match existing names. 
     */
    const string & operator()(unsigned int rowIndex, const string & colName) const
            throw (IndexOutOfBoundsException, NoTableColumnNamesException, TableNameNotFoundException);
    
    /**
     * @name Work on columns.
     *
     * @{
     */

    /**
     * @return The number of columns in this table.
     */
    unsigned int getNumberOfColumns() const { return nCol_; }

    /**
     * @brief Set the column names of this table.
     * 
     * @param colNames The row names.
     * @throw DimensionException If the number of names do not match the number of columns in the table.
     * @throw DuplicatedTableColumnNameException If names are not unique.
     */
    void setColumnNames(const vector<string> & colNames) throw (DimensionException, DuplicatedTableColumnNameException);
    /**
     * @brief Get the column names of this table.
     * 
     * @return The column names of this table.
     * @throw NoTableColumnNamesException If no column names are associated to this table.
     */
    vector<string> getColumnNames() const throw (NoTableColumnNamesException);
    /**
     * @brief Get a given column name.
     * 
     * @param index The index of the column.
     * @return The column name associated to the given column.
     * @throw NoTableColumnNamesException If no column names are associated to this table.
     * @throw IndexOutOfBoundsException If index is >= number of columns.
     */
    string getColumnName(unsigned int index) const throw (NoTableColumnNamesException, IndexOutOfBoundsException);
    
    /**
     * @return true If column names are associated to this table.
     */
    bool hasColumnNames() const { return colNames_!= NULL; }

    /**
     * @return The values in the given column.
     * @param index The index of the column.
     * @throw IndexOutOfBoundsException If index is >= number of columns.
     */
    vector<string> & getColumn(unsigned int index) throw (IndexOutOfBoundsException);
    /**
     * @return The values in the given column.
     * @param index The index of the column.
     * @throw IndexOutOfBoundsException If index is >= number of columns.
     */
    const vector<string> & getColumn(unsigned int index) const throw (IndexOutOfBoundsException);
    
    /**
     * @return The values in the given column.
     * @param colName The name of the column.
     * @throw NoTableColumnNamesException If no column names are associated to this table.
     * @throw TableColumnNameNotFoundException If colName do not match existing column names. 
     */
    vector<string> & getColumn(const string & colName) throw (NoTableColumnNamesException, TableColumnNameNotFoundException);
    /**
     * @return The values in the given column.
     * @param colName The name of the column.
     * @throw NoTableColumnNamesException If no column names are associated to this table.
     * @throw TableColumnNameNotFoundException If colName do not match existing column names. 
     */
    const vector<string> & getColumn(const string & colName) const throw (NoTableColumnNamesException, TableColumnNameNotFoundException);

    /**
     * @brief Tell is a given column exists.
     *
     * @param colName The name of the column to look for.
     * @return true if the column was found, false if not or if there are no column names.
     */
    bool hasColumn(const string & colName) const;

    /**
     * @brief Delete the given column.
     * 
     * @param index The index of the column.
     * @throw IndexOutOfBoundsException If index is >= number of columns.
     */
    void deleteColumn(unsigned int index) throw (IndexOutOfBoundsException);
    
    /**
     * @brief Delete the given column.
     * 
     * @param colName The name of the column.
     * @throw NoTableColumnNamesException If no column names are associated to this table.
     * @throw TableColumnNameNotFoundException If colName do not match existing column names. 
     */
    void deleteColumn(const string & colName) throw (NoTableColumnNamesException, TableColumnNameNotFoundException);
  
    /**
     * @brief Add a new column.
     *
     * @param newColumn The new column values.
     * @throw DimensionException If the number of values does not match the number of rows.
     * @throw TableColumnNamesException If the table has row names.
     */
    void addColumn(const vector<string> & newColumn) throw (DimensionException, TableColumnNamesException);
    /**
     * @brief Add a new column.
     *
     * @param colName   The name of the column.
     * @param newColumn The new column values.
     * @throw DimensionException If the number of values does not match the number of rows.
     * @throw NoTableColumnNamesException If the table does not have row names.
     * @throw DuplicatedTableColumnNameException If colName is already used.
     */
    void addColumn(const string & colName, const vector<string> & newColumn) throw (DimensionException, NoTableColumnNamesException, DuplicatedTableColumnNameException);
    /** @} */
    
    /**
     * @name Work on rows.
     *
     * @{
     */
    
    /**
     * @return The number of rows in this table.
     */
    unsigned int getNumberOfRows() const { return nRow_; }

    /**
     * @brief Set the row names of this table.
     * 
     * @param rowNames The row names.
     * @throw DimensionException If the number of names do not match the number of rows in the table.
     * @throw DuplicatedTableRowNameException If names are not unique.
     */
    void setRowNames(const vector<string> & rowNames) throw (DimensionException, DuplicatedTableRowNameException);

    /**
     * @brief Get the row names of this table.
     * 
     * @return The row names of this table.
     * @throw NoTableRowNamesException If no row names are associated to this table.
     */
    vector<string> getRowNames() const throw (NoTableRowNamesException);

    /**
     * @brief Tell is a given row exists.
     *
     * @param rowName The name of the row to look for.
     * @return true if the row was found, false if not or if there are no row names.
     */
    bool hasRow(const string & rowName) const;

    /**
     * @brief Get a given row name.
     * 
     * @param index The index of the row.
     * @return The row name associated to the given row.
     * @throw NoTableRowNamesException If no row names are associated to this table.
     * @throw IndexOutOfBoundsException If index is >= number of rows.
     */
    string getRowName(unsigned int index) const throw (NoTableRowNamesException, IndexOutOfBoundsException);

    /**
     * @return true If row names are associated to this table.
     */
    bool hasRowNames() const { return rowNames_!= NULL; }

    /**
     * @return A vector which contains a copy  in the given row.
     * @param index The index of the row.
     * @throw IndexOutOfBoundsException If index is >= number of rows.
     */
    vector<string> getRow(unsigned int index) const throw (IndexOutOfBoundsException);

    /**
     * @return A vector which contains a copy  in the given row.
     * @param rowName The name of the row.
     * @throw NoTableRowNamesException If no row names are associated to this table.
     * @throw TableRowNameNotFoundException If rowName do not match existing row names.
     */
    vector<string> getRow(const string & rowName) const throw (NoTableRowNamesException, TableRowNameNotFoundException);

    /**
     * @brief Delete the given row.
     * 
     * @param index The index of the row.
     * @throw IndexOutOfBoundsException If index is >= number of row.
     */
    void deleteRow(unsigned int index) throw (IndexOutOfBoundsException);

    /**
     * @brief Delete the given row.
     * 
     * @param rowName The name of the row.
     * @throw NoTableRowNamesException If no row names are associated to this table.
     * @throw TableRowNameNotFoundException If rowName do not match existing column names. 
     */
    void deleteRow(const string & rowName) throw (NoTableRowNamesException, TableRowNameNotFoundException);
    
    /**
     * @brief Add a new row.
     *
     * @param newRow The new row values.
     * @throw DimensionException If the number of values does not match the number of columns.
     * @throw TableRowNamesException If the table has column names.
     */
    void addRow(const vector<string> & newRow) throw (DimensionException, TableRowNamesException);
    /**
     * @brief Add a new row.
     *
     * @param rowName   The name of the row.
     * @param newRow    The new row values.
     * @throw DimensionException If the number of values does not match the number of columns.
     * @throw NoTableRowNamesException If the table does not have column names.
     * @throw DuplicatedTableRowNameException If rowName is already used.
     */
    void addRow(const string & rowName, const vector<string> & newRow) throw (DimensionException, NoTableRowNamesException, DuplicatedTableRowNameException);
    /** @} */

  public:

    /**
     * @brief Read a table form a stream in CSV-like format.
     *
     * The number of rows is given by the second line in the file.
     * By default, if the first line as one column less than the second one,
     * the first line is taken as column names, and the first column as row names.
     * Otherwise, no column names and no row names are specified, unless
     * explicitely precised by the user.
     * 
     * @param in       The input stream.
     * @param sep      The column delimiter.
     * @param header   Tell if the first line must be used as column names, otherwise use default.
     * @param rowNames Use a column as rowNames. If positive, use the specified column to compute rownames, otherwise use default;
     * @return         A pointer toward a new DataTable object.
     */
    static DataTable * read(istream & in, const string & sep = "\t", bool header=true, int rowNames=-1)
      throw (DimensionException, IndexOutOfBoundsException, DuplicatedTableRowNameException);

    /**
     * @brief Write a DataTable object to stream in CVS-like format.
     * 
     * @param data     The table to write.
     * @param out      The output stream.
     * @param sep      The column delimiter.
     */
    static void write(const DataTable & data, ostream & out, const string & sep = "\t");
};

} //end of namespace bpp.

#endif //_DataTable_H_

