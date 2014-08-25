#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <queue>
#include <algorithm>

using namespace std;

set<int> all_possible_values = {1,2,3,4,5,6,7,8,9};

class sudoku_cell {

  set<int> possible_values;
  int cell_value = 0;
  std::pair<size_t,size_t> coors;
public:
  void print_coors() { cout <<" "<<coors.first<<","<<coors.second<<" "; }
  const pair<size_t,size_t>& return_coors() const { return coors; } 
  bool value_included_in_possible_values_set( int value ) const;
  const set<int>& return_possible_values_set() const; 
  int count_remaining_possible_values() const; 
  void exclude_value( int value );
  void set_value( int value );
  void print_cell_value() const; 
  const int& return_cell_value() const; 
  sudoku_cell( int input, size_t i, size_t j );
  ~sudoku_cell() { ; }
}; 


bool sudoku_cell::value_included_in_possible_values_set( int value ) const {
  
  for( int val : possible_values )
    if( val == value ) return true;    
  return false;  
}  

void sudoku_cell::set_value( int value ) {

  cell_value = value;
  possible_values.clear();
}  

const set<int>& sudoku_cell::return_possible_values_set() const { 

  return possible_values;  
}  

void sudoku_cell::exclude_value( int value ) {

  if( cell_value != 0 ) return;

  auto res = possible_values.erase(value); 
  if( res == 1 && possible_values.size() == 1 ) 
    set_value( *(possible_values.cbegin()) );   
}  

sudoku_cell::sudoku_cell( int input, size_t i, size_t j ) { 
 
  if( input != -1 ) cell_value = input; 
  else possible_values = all_possible_values;
  this-> coors = std::pair<size_t,size_t>(i,j);
  //this-> j = j;
}

const int& sudoku_cell::return_cell_value() const {
  // returned values are 0 for undefined valued cells, or the cell's final value otherwise
  return cell_value;
}  

void sudoku_cell::print_cell_value() const {

  if(possible_values.size() > 1) { 
    
    for( int val : possible_values ) cout<<val;
    for(size_t i=0; i<9-possible_values.size(); i++) cout <<" ";

  } else cout << cell_value<<"        "; //only one value to print
}  

int sudoku_cell::count_remaining_possible_values() const { 

  return possible_values.size();
}  


// This class is used to sort the various sudoku cells. The metric used here is the
// number of possible values that the cell can have. The less possible values a cell
// has, means that it probably will require less work in order to find it's final value.
// That's why it should go higher in the queue, while the cells that require more work are
// left at lower levels of the container of the priority queue.
class unsolved_sudoku_cell_compare_method {

public:
  bool operator ()(sudoku_cell* const & first_cell, sudoku_cell* const & second_cell) const {

    return first_cell-> count_remaining_possible_values() < // ascending order sorting 
           second_cell-> count_remaining_possible_values();  
  }  
};  


class sudoku_solver {

  struct candidate_value_struct{
     int value=0;
     set<pair<size_t,size_t>> appearing_cells={};     
  }; 
  size_t total_dangling_values = 9*9*9;
  // The 3 arguments are: object type, container of object, compare method
  priority_queue<sudoku_cell*, vector<sudoku_cell*>, unsolved_sudoku_cell_compare_method> unsolved_cells_queue;
  set<int> gathered_values = {};
  int cell_inputs[9][9];
  sudoku_cell* sudoku_array[9][9];
  void observe_cell_row( sudoku_cell* const & examined_cell );
  void observe_cell_column( sudoku_cell* const & examined_cell );
  void observe_cell_subblock( sudoku_cell* const & examined_cell );
  void try_further_possible_values_eliminaton( sudoku_cell* const & examined_cell );  
  tuple<bool,int,int,pair<size_t,size_t>> 
  seeked_possible_value_is_in_single_row_or_column_of_subblock( const int & seeked_value, 
                                                                const size_t & subblock_row, 
                                                                const size_t & subblock_column );
  void examine_for_pairs( sudoku_cell* const & examined_cell );
  void examine_for_triplets( sudoku_cell* const & examined_cell );
  void search_for_easy_naked_pairs( const vector<sudoku_cell*> & collected_cells );
  void search_for_naked_pairs( vector<struct candidate_value_struct> & collected_candidate_value_structs,
                               const set<pair<size_t,size_t>> & collected_coors );
  void search_for_easy_naked_triplets( const vector<sudoku_cell*> & collected_cells );
  void search_for_hidden_triplets( vector<struct candidate_value_struct> & collected_candidate_value_structs,
                                   const set<pair<size_t,size_t>> & collected_coors );
  void initialize_collected_candidate_value_structs( vector<struct candidate_value_struct> & data_vec, 
                                                     set<pair<size_t,size_t>> & collected_coors, 
                                                     vector<sudoku_cell*> & collected_cells ); 
  size_t count_dangling_values();
public:  
  sudoku_solver() { // Initialize the queue in the constructor 
    unsolved_cells_queue = move(priority_queue<sudoku_cell*, vector<sudoku_cell*>, unsolved_sudoku_cell_compare_method>()); 
  }
  void print_grid();
  void solve();
  void collect_input( char* input_file );
  ~sudoku_solver() {         
    for(size_t i=0; i<9; i++)
      for(size_t j=0; j<9; j++)
        if(sudoku_array[i][j] != nullptr) sudoku_array[i][j]->~sudoku_cell(); 
  }
};  

size_t sudoku_solver::count_dangling_values() {

  size_t total = 0; 
  for(size_t i=0; i<9; i++)
    for(size_t j=0; j<9; j++)
      total += sudoku_array[i][j]-> return_possible_values_set().size();
  return total;
}  

void sudoku_solver::print_grid() { 

  for(size_t i=0; i<9; i++){
     for(size_t j=0; j<9; j++) {
        sudoku_array[i][j]-> print_cell_value(); cout<<" ";
     } cout << endl;
  }   
}

void sudoku_solver::observe_cell_row( sudoku_cell* const & examined_cell ) {

  if( examined_cell-> return_cell_value() != 0 ) return; // if final cell value is set, return
  gathered_values.clear();

  auto row_numb = (examined_cell-> return_coors()).first; 
  auto column_numb = (examined_cell-> return_coors()).second; 

  for(size_t i=0; i<9; i++) // gather all the possible values that are yet to appear on this row 
    if( i!=column_numb ) {
      auto possible_cell_values_set = sudoku_array[row_numb][i]-> return_possible_values_set(); 
      for( int val : possible_cell_values_set ) gathered_values.insert( val ); 
    }

  cout <<"Values that haven't appeared yet in row "<< row_numb 
       << ",excluding cell ["<<row_numb<<"]["<<column_numb<<"] are:"; 
  for( int val : gathered_values ) cout << val <<","; cout <<endl;

  for( int possible_val : examined_cell-> return_possible_values_set() )
    // if possible_val is not found anywhere else in this row, then 
    if( gathered_values.end() == std::find(gathered_values.begin(), gathered_values.end(), possible_val) ) {
       examined_cell-> set_value( possible_val ); break; 
    }   
}
 

void sudoku_solver::observe_cell_column( sudoku_cell* const & examined_cell ) {

  if( examined_cell-> return_cell_value() != 0 ) return; // Redundant here
  gathered_values.clear();

  auto row_numb = (examined_cell-> return_coors()).first; 
  auto column_numb = (examined_cell-> return_coors()).second; 

  for(size_t i=0; i<9; i++) // gather all the possible values that are yet to appear on this row 
      if( i!=row_numb ) {
        auto possible_cell_values_set = sudoku_array[i][column_numb]-> return_possible_values_set(); 
        for( int val : possible_cell_values_set ) gathered_values.insert( val ); 
      }

  cout <<"Values that haven't appeared yet in column "<< column_numb 
       << ",excluding cell ["<<row_numb<<"]["<<column_numb<<"] are:"; 
  for( int val : gathered_values ) cout << val <<","; cout <<endl;
  
  for( int possible_val : examined_cell-> return_possible_values_set() )
    // if possible_val is not found anywhere else in this row, then 
    if( gathered_values.end() == std::find(gathered_values.begin(), gathered_values.end(), possible_val) ) {
       examined_cell-> set_value( possible_val ); break; 
    } 
}  


void sudoku_solver::observe_cell_subblock( sudoku_cell* const & examined_cell ) {

  if( examined_cell-> return_cell_value() != 0 ) return;
  gathered_values.clear();

  size_t row_numb = examined_cell->return_coors().first;   
  size_t column_numb = examined_cell->return_coors().second;
  size_t subblock_row = (row_numb/3)*3, subblock_column = (column_numb/3)*3;

  for(size_t i = subblock_row; i<subblock_row+3; i++) 
    for(size_t j = subblock_column; j<subblock_column+3; j++)   
      if( i == row_numb && j == column_numb ) ;   
      else { 
        auto possible_cell_values_set = sudoku_array[i][j]-> return_possible_values_set(); 
        for( int val : possible_cell_values_set ) gathered_values.insert( val ); 
      }

  cout <<"Values that haven't appeared yet in subblock "<< subblock_row 
       << ","<<subblock_column<< ",excluding cell ["<<row_numb<<"]["<<column_numb<<"] are:"; 
  for( int val : gathered_values ) cout << val <<","; cout <<endl;
 
  for( int possible_val : examined_cell-> return_possible_values_set() )
    // if possible_val is not found anywhere else in this row, then 
    if( gathered_values.end() == std::find(gathered_values.begin(), gathered_values.end(), possible_val) ) {
       examined_cell-> set_value( possible_val ); break; 
    }  
}  


// For the cells without a final value yet, query and see which ones have seeked_value in their possible_values set. 
// From examining the coordinates, we can see if the seeked_value is expected to be on a single row or column, or
// if it's expected to appear in more than one rows and/or columns
tuple<bool,int,int,pair<size_t,size_t>> 
sudoku_solver::seeked_possible_value_is_in_single_row_or_column_of_subblock(const int & seeked_value, 
                                                                            const size_t & subblock_row, 
                                                                            const size_t & subblock_column ) {

  vector<pair<size_t,size_t>> coordinates_containing_seeked_value = {};
  int row_of_interest = -1; int column_of_interest = -1;

  for(size_t i = subblock_row; i<subblock_row+3; i++) 
    for(size_t j = subblock_column; j<subblock_column+3; j++) {
 
      if( sudoku_array[i][j]-> return_cell_value() == seeked_value )
        // This subblock has the seeked_value already present as a permanent value in cell [i][j],
        // so no need to continue 
        return tuple<bool,int,int,pair<size_t,size_t>>(false,-1,-1,pair<size_t,size_t>(0,0)); 
            
      if( sudoku_array[i][j]-> return_cell_value() == 0 ) // cells with no permanent value yet 
        if( sudoku_array[i][j]-> value_included_in_possible_values_set( seeked_value ) ) // if condition is true
               coordinates_containing_seeked_value.push_back( pair<size_t,size_t>(i,j) ); 
    }

  // Comment this case more!! 
  // In this subblock, the value seeked_value will be found in ONLY one place. However, that value hasn't been
  // assigned yet to sudoku_cell::cell_value, but exists in the sudoku_cell::possible_values set.
  // It is tempting to set that cell's value here. HOWEVER, we don't, to be more clear in the code
  // and avoid confusions.
  if( coordinates_containing_seeked_value.size() == 1 )
    return tuple<bool,int,int,pair<size_t,size_t>>(false,-1,-1,pair<size_t,size_t>(0,0)); 

  if( coordinates_containing_seeked_value.size() == 0 ) { 
    cout<<"In "<< subblock_row <<","<<subblock_column<<" ";
    cout<<"Cannot find "<<seeked_value<<" anywhere!!"; exit(0);     
  }

  // In this subblock, we've retrieved the cells that could contain the value seeked_value.
  // By examining the coordinates, we can see if a value is expected to reside on a single row or column,
  // or maybe in more than one row and/or column. If the value is expected to reside in a single row ONLY, 
  // or a single column ONLY, that information is sent back.
  // WE now have all the coordinates. 

  // We care to see if all coordinates have a common x or y parameter. That would mean that the seeked_value lies
  // either on row x or on column y.
  // To achieve that, we sort all (x,y) first using x, and then using y as parameter of interest.
  // If after the sort, the first sorted element AND the last sorted element have the same x or y parameter, 
  // then seeked_value definitely will reside in row x or columny in this subblock.

  std::sort( coordinates_containing_seeked_value.begin(), coordinates_containing_seeked_value.end(), 
             [](const pair<size_t,size_t>& coor1, const pair<size_t,size_t>& coor2 )-> bool {return coor1.first<coor2.first;});

  if( (*coordinates_containing_seeked_value.cbegin()).first == 
      (*coordinates_containing_seeked_value.crbegin()).first ) {
 
      cout << "In subblock ["<<subblock_row<<"]["<<subblock_column<<"], value "
           <<seeked_value<<" can be found ONLY in row "
           <<(*coordinates_containing_seeked_value.crbegin()).first<<endl;

      row_of_interest = (*coordinates_containing_seeked_value.crbegin()).first;
  }  

  // sort (x,y) coordinates using y as parameter of interest. 
  std::sort( coordinates_containing_seeked_value.begin(), coordinates_containing_seeked_value.end(), 
             []( const pair<size_t,size_t>& coor1, const pair<size_t,size_t>& coor2 )-> bool {return coor1.second<coor2.second;});

  if( (*coordinates_containing_seeked_value.cbegin()).second == 
      (*coordinates_containing_seeked_value.crbegin()).second ) {
 
      cout << "In subblock ["<<subblock_row<<"]["<<subblock_column<<"], value "
           <<seeked_value<<" can be found ONLY in column "
           <<(*coordinates_containing_seeked_value.crbegin()).second<<endl;

      column_of_interest = (*coordinates_containing_seeked_value.crbegin()).second;
  }

  if( row_of_interest != -1 && column_of_interest != -1 ) { cout << "Problem with seeked_value!"; exit(0); }

  // return value is false, if we haven't discovered a row or column of interest. 
  return tuple<bool,int,int,pair<size_t,size_t>>( row_of_interest != -1 || column_of_interest != -1,
                                                  row_of_interest, column_of_interest,
                                                  pair<size_t,size_t>(subblock_row, subblock_column));
}


void sudoku_solver::try_further_possible_values_eliminaton( sudoku_cell* const & examined_cell ) {

  if( examined_cell-> return_cell_value() != 0 ) return;
  gathered_values.clear();

  size_t cell_row = examined_cell->return_coors().first;   
  size_t cell_column = examined_cell->return_coors().second;

  set<int> my_set = examined_cell-> return_possible_values_set();
  for( auto possible_value : my_set )
  for( size_t row=0; row<9; row=row+3 )
    for( size_t column=0; column<9; column=column+3 ) {

      if( row == (cell_row/3)*3 && column == (cell_column/3)*3 ) ; // excluding examined_cell's subblock
      else { 
        /*return type is tuple<bool,int,int,pair<size_t,size_t>>*/
        auto subblock_result = seeked_possible_value_is_in_single_row_or_column_of_subblock( possible_value, 
                                                                                             row, column ); 
        if(std::get<0>(subblock_result)==true) { // if this is false, then we don't have enough useful info to proceed 
          if( std::get<1>(subblock_result) == cell_row ) { // either cell_row or cell_column will be equal
         
            cout <<"Removing possible value "<<possible_value<<" from cell ["
                 <<cell_row<<"]["<<cell_column<<"], it is to be expected in row "<<cell_row
                 <<" in subblock ["<<std::get<3>(subblock_result).first
                 <<"]["<< std::get<3>(subblock_result).second<<"]"<<endl;

            examined_cell-> exclude_value( possible_value );

          } else if( std::get<2>(subblock_result) == cell_column ) { 
          
            cout <<"Removing possible value "<<possible_value<<" from cell ["
                 <<cell_row<<"]["<<cell_column<<"], it is to be expected in column "<<cell_column
                 <<" in subblock ["<<std::get<3>(subblock_result).first
                 <<"]["<< std::get<3>(subblock_result).second<<"]"<<endl;

            examined_cell-> exclude_value( possible_value );
          }
        }  
      }     
    }
}



void sudoku_solver::initialize_collected_candidate_value_structs( vector<struct candidate_value_struct> & data_vec,
                                                                  set<pair<size_t,size_t>> & collected_coors,
                                                                  vector<sudoku_cell*> & collected_cells ) {
  collected_cells.clear(); 
  collected_coors.clear();
  data_vec.clear(); // clear previous values, if any
  for(size_t i=0; i<9; i++) { // Initialization 
    struct candidate_value_struct i_value_struct;  
    i_value_struct.value = i+1;
    data_vec.push_back( std::move(i_value_struct) ); 
  } 
}  


void sudoku_solver::examine_for_pairs( sudoku_cell* const & examined_cell ) {

  // Index is used for representing the value-> index 3 corresponds to value 4
  vector<struct candidate_value_struct> collected_candidate_value_structs = {};
  set<std::pair<size_t,size_t>> collected_coors = {};
  vector<sudoku_cell*> collected_cells = {};

  size_t cell_row = examined_cell->return_coors().first;   
  size_t cell_column = examined_cell->return_coors().second;

  initialize_collected_candidate_value_structs( collected_candidate_value_structs, collected_coors, collected_cells );

  for(size_t j=0; j<9; j++) collected_cells.push_back( sudoku_array[cell_row][j] );
  search_for_easy_naked_pairs( collected_cells );

  for(size_t j=0; j<9; j++) { // examining row cell_row  
 
    collected_coors.insert( pair<size_t,size_t>(cell_row,j) ); 
    int cell_value = sudoku_array[cell_row][j]-> return_cell_value();
    if( cell_value != 0 ) { // this cell has a final value 
        collected_candidate_value_structs[ cell_value-1 ].appearing_cells = { pair<size_t,size_t>(cell_row,j) };    
    } else {
      for( int possible_value : sudoku_array[cell_row][j]-> return_possible_values_set() ) 
          collected_candidate_value_structs[ possible_value-1 ].appearing_cells.insert(pair<size_t,size_t>(cell_row,j));     
    }
  } // data collection from row has ended.  
  search_for_naked_pairs( collected_candidate_value_structs, collected_coors ); // search for naked pairs in row cell_row
  // processing data of row has ended 

  initialize_collected_candidate_value_structs( collected_candidate_value_structs, collected_coors, collected_cells );

  for(size_t i=0; i<9; i++) collected_cells.push_back( sudoku_array[i][cell_column] );
  search_for_easy_naked_pairs( collected_cells );

  for(size_t i=0; i<9; i++) { // data collection from column
    
    collected_coors.insert( pair<size_t,size_t>(i,cell_column) ); 
    int cell_value = sudoku_array[i][cell_column]-> return_cell_value();
    if( cell_value != 0 ) { // this cell has a final value 
        collected_candidate_value_structs[ cell_value-1 ].appearing_cells = { pair<size_t,size_t>(i,cell_column) };    
    } else {
      for( int possible_value : sudoku_array[i][cell_column]-> return_possible_values_set() ) 
          collected_candidate_value_structs[ possible_value-1 ].appearing_cells.insert(pair<size_t,size_t>(i,cell_column));
    }
  } 
  search_for_naked_pairs( collected_candidate_value_structs, collected_coors ); 
  // search for naked pairs in column cell_column has ended 

  initialize_collected_candidate_value_structs( collected_candidate_value_structs, collected_coors, collected_cells ); 

  size_t subblock_row = (cell_row/3)*3, subblock_column = (cell_column/3)*3;

  for(size_t i = subblock_row; i<subblock_row+3; i++) 
    for(size_t j = subblock_column; j<subblock_column+3; j++) collected_cells.push_back( sudoku_array[i][j] );
 
  search_for_easy_naked_pairs( collected_cells );

  for(size_t i = subblock_row; i<subblock_row+3; i++) 
    for(size_t j = subblock_column; j<subblock_column+3; j++) { 
 
      collected_coors.insert( pair<size_t,size_t>(i,j) ); 
      int cell_value = sudoku_array[i][j]-> return_cell_value();
      if( cell_value != 0 ) { // this cell has a final value 
        collected_candidate_value_structs[ cell_value-1 ].appearing_cells = { pair<size_t,size_t>(i,j) };    
      } else {
         for( int possible_value : sudoku_array[i][j]-> return_possible_values_set() ) 
           collected_candidate_value_structs[ possible_value-1 ].appearing_cells.insert(pair<size_t,size_t>(i,j));  
      }
    } // data collection from subblock  
  search_for_naked_pairs( collected_candidate_value_structs, collected_coors );
}  

void sudoku_solver::search_for_easy_naked_pairs( const vector<sudoku_cell*> & collected_cells ) {

  // This function identifies cells that have IDENTICAL possible_values sets, and which contain 2 elements only
  // If so, these two elements will be removed from the possible_values sets of the rest
  // sudoku sells, contained inside the vector collected_cells 

  // Now, we compare the possible_value sets A,B,...N. IF A==B, or ... A==N, and have a size of 2, then
  // they are a naked pair. We try all combos  
  for(size_t i=0; i<collected_cells.size(); i++)
    for(size_t j=i+1; j<collected_cells.size(); j++) // lambda function
     [&]( const set<int> & possible_value_set_A, const set<int> & possible_value_set_B )-> void {  
 
        if( possible_value_set_A.size() == 2 && possible_value_set_A == possible_value_set_B ) { 

          // the sets are identical (remember they are stored in an RB tree), the == operator works on sets  
          cout <<"***Easy naked pair! Values ";
          for( int possible_value : possible_value_set_A ) cout << possible_value<<",";
          cout <<" will appear as a pair in "; collected_cells[i]-> print_coors();
          cout <<" and "; collected_cells[j]-> print_coors(); 
          cout << endl;              

          // remove these two possible values from the rest of the sudoku cells     
          for( sudoku_cell* const ptr : collected_cells ) 
            if( ptr != collected_cells[i] && ptr != collected_cells[j] ) 
              for( int possible_value : possible_value_set_A ) ptr-> exclude_value( possible_value ); 
        }

     } ( collected_cells[i]->return_possible_values_set(), collected_cells[j]->return_possible_values_set() );

}


// Using this function, in the data set that we have now, we try to find hidden pairs 
void sudoku_solver::search_for_naked_pairs( vector<struct candidate_value_struct> & collected_candidate_value_structs,
                                            const set<pair<size_t,size_t>> & collected_coors ) {

  // Finds all possible values that appear exactly 2 times. The first element is found at the beginning, 
  // and the last element is pointed one position before what's pointed by end_iter 
  auto end_iter = std::remove_if( collected_candidate_value_structs.begin(),
                                  collected_candidate_value_structs.end(), 
                                  // If this function returns true, the element is removed 
                                  []( struct candidate_value_struct & cvs )-> bool { return cvs.appearing_cells.size()!=2; });

  // Chopping the rest elements, keeping only those with numbor of appearances==2
  collected_candidate_value_structs.erase( end_iter, collected_candidate_value_structs.end() );

  if( collected_candidate_value_structs.size() <= 1 ) return;

  cout << collected_candidate_value_structs.size() <<" possible values can appear 2 times here"<<endl;
  // A naked pair MUST have only two cells that can appear. So, each element of the naked pair will have 2 cells
  // where it can potentially be. IF these two cells are the same, then we know the 2 cells that ONLY there can
  // that naked pair reside
  for( auto cvs : collected_candidate_value_structs ) {
    cout<<"Value "<< cvs.value << " appears in: ";
    for( auto cells : cvs.appearing_cells ) cout<<"["<<cells.first<<"]["<<cells.second<<"], ";
  } cout<< endl;

  // We now have 2+ structs A,B,C...,N and we must compare each pair: (A,B),(A,C),...(A,N),(B,C),...(B,N),...(N-1,N)
  // Try all pair combinations of A,B,C,..,N -> (A,B),(A,C),...,(A,N),(B,C),...,(B,N),...,(N-1,N)
  // Remember that A,B,..,N are objects of type struct candidate_value_struct
  for(size_t i=0; i<collected_candidate_value_structs.size(); i++)
    for(size_t j=i+1; j<collected_candidate_value_structs.size(); j++) // lambda function
     [&]( set<pair<size_t,size_t>> & coor_set_A, set<pair<size_t,size_t>> & coor_set_B )-> void {  

       // We examine the coordinates that value_1 and value_2 appear. Remember both will have 2 coordinates only!  
       // Say value_1 has coordinates (x1,y1),(z1,l1), and value_2 has (x2,y2),(z2,l2), remember that both are sorted
       // In naked pairs, BOTH coordinates must be the same, otherwise we don't have a naked pair 
       if( coor_set_A == coor_set_B ) {   

         cout <<"***Naked pair! Values "<<collected_candidate_value_structs[i].value<<" and "
              <<collected_candidate_value_structs[j].value<<" will appear as a pair in :";
              for( auto & coor : coor_set_A ) cout<<"["<<coor.first<<"]["<<coor.second<<"], "; cout<<endl;

         for( auto & coor : coor_set_A )
         for( int k=1; k<10; k++ ) 
           if( k != collected_candidate_value_structs[i].value && k != collected_candidate_value_structs[j].value ) 
             sudoku_array[coor.first][coor.second]-> exclude_value(k);
  
         //Now, must also remove the possible_values pair from the other examined cells
         set<pair<size_t,size_t>> diff_set = collected_coors; // constructing the difference set 
         for( auto coor : coor_set_A ) diff_set.erase(coor); 

         for( auto & coll_coor : diff_set ) { 
             sudoku_array[coll_coor.first][coll_coor.second]-> exclude_value( collected_candidate_value_structs[i].value );
             sudoku_array[coll_coor.first][coll_coor.second]-> exclude_value( collected_candidate_value_structs[j].value );
         } 
       }
     } ( collected_candidate_value_structs[i].appearing_cells, collected_candidate_value_structs[j].appearing_cells );  
}  


void sudoku_solver::examine_for_triplets( sudoku_cell* const & examined_cell ) {

// Index is used for representing the value-> index 3 corresponds to value 4
  vector<struct candidate_value_struct> collected_candidate_value_structs = {};
  set<std::pair<size_t,size_t>> collected_coors = {};
  vector<sudoku_cell*> collected_cells = {};

  size_t cell_row = examined_cell->return_coors().first;   
  size_t cell_column = examined_cell->return_coors().second;

  initialize_collected_candidate_value_structs( collected_candidate_value_structs, collected_coors, collected_cells );

  for(size_t j=0; j<9; j++) collected_cells.push_back( sudoku_array[cell_row][j] );
  search_for_easy_naked_triplets( collected_cells );

  for(size_t j=0; j<9; j++) { // examining row cell_row  
 
    collected_coors.insert( pair<size_t,size_t>(cell_row,j) ); 
    int cell_value = sudoku_array[cell_row][j]-> return_cell_value();
    if( cell_value != 0 ) { // this cell has a final value 
        collected_candidate_value_structs[ cell_value-1 ].appearing_cells = { pair<size_t,size_t>(cell_row,j) };    
    } else {
      for( int possible_value : sudoku_array[cell_row][j]-> return_possible_values_set() ) 
          collected_candidate_value_structs[ possible_value-1 ].appearing_cells.insert(pair<size_t,size_t>(cell_row,j)); 
    }
  } // data collection from row has ended.  
  search_for_hidden_triplets( collected_candidate_value_structs, collected_coors );
  // processing data of row has ended 

  initialize_collected_candidate_value_structs( collected_candidate_value_structs, collected_coors, collected_cells );

  for(size_t i=0; i<9; i++) collected_cells.push_back( sudoku_array[i][cell_column] );
  search_for_easy_naked_triplets( collected_cells );

  for(size_t i=0; i<9; i++) { // data collection from column
    
    collected_coors.insert( pair<size_t,size_t>(i,cell_column) ); 
    int cell_value = sudoku_array[i][cell_column]-> return_cell_value();
    if( cell_value != 0 ) { // this cell has a final value 
        collected_candidate_value_structs[ cell_value-1 ].appearing_cells = { pair<size_t,size_t>(i,cell_column) };    
    } else {
      for( int possible_value : sudoku_array[i][cell_column]-> return_possible_values_set() ) 
          collected_candidate_value_structs[ possible_value-1 ].appearing_cells.insert(pair<size_t,size_t>(i,cell_column));    
    }
  } 
  search_for_hidden_triplets( collected_candidate_value_structs, collected_coors ); 
  // search for triplets in column cell_column has ended 

  initialize_collected_candidate_value_structs( collected_candidate_value_structs, collected_coors, collected_cells ); 

  size_t subblock_row = (cell_row/3)*3, subblock_column = (cell_column/3)*3;

  for(size_t i = subblock_row; i<subblock_row+3; i++) 
    for(size_t j = subblock_column; j<subblock_column+3; j++) collected_cells.push_back( sudoku_array[i][j] );
 
  search_for_easy_naked_triplets( collected_cells );

  for(size_t i = subblock_row; i<subblock_row+3; i++) 
    for(size_t j = subblock_column; j<subblock_column+3; j++) { 
 
      collected_coors.insert( pair<size_t,size_t>(i,j) ); 
      int cell_value = sudoku_array[i][j]-> return_cell_value();
      if( cell_value != 0 ) { // this cell has a final value 
        collected_candidate_value_structs[ cell_value-1 ].appearing_cells = { pair<size_t,size_t>(i,j) };    
      } else {
          for( int possible_value : sudoku_array[i][j]-> return_possible_values_set() ) 
            collected_candidate_value_structs[ possible_value-1 ].appearing_cells.insert(pair<size_t,size_t>(i,j));
      }
    } // data collection from subblock  
  search_for_hidden_triplets( collected_candidate_value_structs, collected_coors );  
}


void sudoku_solver::search_for_easy_naked_triplets( const vector<sudoku_cell*> & collected_cells ) {
  
  // This function identifies cells that have IDENTICAL possible_values sets, and which contain 2 or 3 elements only
  // If so, these two elements will be removed from the possible_values sets of the rest
  // sudoku sells, contained inside the vector collected_cells 

  // Now, we compare the possible_value sets A,B,...N. IF A==B, or ... A==N, and have a size of 2, then
  // they are a naked pair. We try all combos  
  for(size_t i=0; i<collected_cells.size(); i++)
    for(size_t j=i+1; j<collected_cells.size(); j++) // lambda function
      for(size_t k=j+1; k<collected_cells.size(); k++) // lambda function
        [&]( const set<int> & possible_vals_A, 
             const set<int> & possible_vals_B,
             const set<int> & possible_vals_C )-> void {  

          if( possible_vals_A.size() == 3 && possible_vals_B.size() == 3 && possible_vals_C.size() == 3 ) 
          if( possible_vals_A == possible_vals_B && possible_vals_B == possible_vals_C ) { 

            // the sets are identical (remember they are stored in an RB tree), the == operator works on sets  
            cout <<"***Easy naked triplet! Values ";
            for( int possible_value : possible_vals_A ) cout << possible_value<<",";
            cout <<" will appear as a triplet in "; collected_cells[i]-> print_coors();
            cout <<" and "; collected_cells[j]-> print_coors(); 
            cout <<" and "; collected_cells[k]-> print_coors();
            cout << endl;

            // remove these three possible values from the rest of the sudoku cells     
            for( sudoku_cell* const ptr : collected_cells ) 
              if( ptr != collected_cells[i] && ptr != collected_cells[j] && ptr != collected_cells[k] ) 
                for( int possible_value : possible_vals_A ) ptr-> exclude_value( possible_value ); 
          }

        } ( collected_cells[i]->return_possible_values_set(), 
            collected_cells[j]->return_possible_values_set(),
            collected_cells[k]->return_possible_values_set() );
}


void sudoku_solver::search_for_hidden_triplets( vector<struct candidate_value_struct> & collected_candidate_value_structs,
                                                const set<pair<size_t,size_t>> & collected_coors ) {

  // Finds all possible values that appear exactly 2 or 3 times. The first element is found at the beginning, 
  // and the last element is pointed one position before what's pointed by end_iter 
  auto end_iter = std::remove_if( collected_candidate_value_structs.begin(),
                                  collected_candidate_value_structs.end(), 
                                  // If this function returns true, the element is removed 
                                  [&]( const struct candidate_value_struct & cvs )-> bool { 
                                    return cvs.appearing_cells.size()==1 || cvs.appearing_cells.size()>3; 
                                  }); 
  // Chopping the rest elements, keeping only those with number of appearances==2 or 3
  collected_candidate_value_structs.erase( end_iter, collected_candidate_value_structs.end() );

  if( collected_candidate_value_structs.size() <= 2 ) return;

  cout << collected_candidate_value_structs.size() <<" possible values can appear 2 or 3 times here :"; 
  for( auto val : collected_candidate_value_structs ) cout << val.value<<", "; cout<<endl;

  // We are to compare a triplet of coordinate sets, set_A, set_B, set_C, all of which have
  // 2 to 3 elements. We compute the union set : set_A U set_B U set_C, which must have 
  // 3 elements exactly, if it's 3 coordinates are to be a hidden triplet 
  for(size_t i=0; i<collected_candidate_value_structs.size(); i++)
    for(size_t j=i+1; j<collected_candidate_value_structs.size(); j++) 
      for(size_t k=j+1; k<collected_candidate_value_structs.size(); k++) 
        [&]( const set<pair<size_t,size_t>> & coor_set_A,
             const set<pair<size_t,size_t>> & coor_set_B,
             const set<pair<size_t,size_t>> & coor_set_C )-> void {  

            set<pair<size_t,size_t>> coors_union_result = {}; 
            // Throw everything in the set, duplicates are not copied, so this way we get a union
            for( auto val : coor_set_A ) coors_union_result.insert(val);
            for( auto val : coor_set_B ) coors_union_result.insert(val);
            for( auto val : coor_set_C ) coors_union_result.insert(val);
 
            if( coors_union_result.size() != 3 ) return; // The union must be exaclty 3 coordinates! 
 
            cout <<"***Hidden Triplet! Values "<<collected_candidate_value_structs[i].value<<", "
                 <<collected_candidate_value_structs[j].value<<", "<< collected_candidate_value_structs[k].value
                 <<" will appear as triplets in "; for( pair<size_t,size_t> coor : coors_union_result ) 
                                                     cout<<"["<<coor.first<<"]["<<coor.second<<"] ";cout<<endl; 
                              
            for( pair<size_t,size_t> coor : coors_union_result )  
              for( int l=1; l<10; l++ ) // remove all but the 3 found values.
                if( l != collected_candidate_value_structs[i].value
                     && l != collected_candidate_value_structs[j].value      
                      && l != collected_candidate_value_structs[k].value ) 
                    sudoku_array[coor.first][coor.second]-> exclude_value( l );   
                 
            //Now, must also remove the possible_vals_set from the remaining examined cells                
            set<pair<size_t,size_t>> rest_coors = collected_coors;
            for( pair<size_t,size_t> excluded_coor : coors_union_result ) rest_coors.erase( excluded_coor );

            for( pair<size_t,size_t> coor : rest_coors ) { 
                  sudoku_array[coor.first][coor.second]-> exclude_value( collected_candidate_value_structs[i].value );
                  sudoku_array[coor.first][coor.second]-> exclude_value( collected_candidate_value_structs[j].value );
                  sudoku_array[coor.first][coor.second]-> exclude_value( collected_candidate_value_structs[k].value ); 
            } cout<<endl;

         }( collected_candidate_value_structs[i].appearing_cells, 
            collected_candidate_value_structs[j].appearing_cells,
            collected_candidate_value_structs[k].appearing_cells ); 
}


void sudoku_solver::solve() {

  auto examined_unsolved_cells_queue = priority_queue<sudoku_cell*, vector<sudoku_cell*>, 
                                                      unsolved_sudoku_cell_compare_method>();
                                        
  while( unsolved_cells_queue.size() > 0 ) {  

    while( unsolved_cells_queue.size() > 0 ) {

     total_dangling_values = count_dangling_values();

     for(size_t i=0; i<9; i++)
     for(size_t j=0; j<9; j++) 
      if( sudoku_array[i][j]-> return_cell_value() != 0 ) {
        [&] ( int final_cell_value ) { // examining row
          // Excluding final_cell_value from ALL the possible_values_set of the other cells in the row
          for(size_t k=0; k<9; k++) sudoku_array[i][k]-> exclude_value( final_cell_value );
        } ( sudoku_array[i][j]-> return_cell_value() );
        
        [&] ( int final_cell_value ) { //examining column

          for(size_t l=0; l<9; l++) sudoku_array[l][j]-> exclude_value( final_cell_value ); 
        } ( sudoku_array[i][j]-> return_cell_value() );

        [&] ( int final_cell_value ) { //examining subblock

          size_t subblock_row = (i/3)*3, subblock_column = (j/3)*3;
          for(size_t k = subblock_row; k<subblock_row+3; k++) 
            for(size_t l = subblock_column; l<subblock_column+3; l++) 
              sudoku_array[k][l]-> exclude_value( final_cell_value ); 
        } ( sudoku_array[i][j]-> return_cell_value() );
      }
  
      auto examined_cell = unsolved_cells_queue.top();       

      // If a method fails to solve the sudoku further, we try the next one
      // If a method succeeds, then the boolean condition won't hold and the remaining 
      // solving methods won't be used  
      if( total_dangling_values == count_dangling_values() ) {

        cout<<"processing cell "; examined_cell-> print_coors(); cout<<endl;
        observe_cell_row( examined_cell ); 
        observe_cell_column( examined_cell );  
        observe_cell_subblock( examined_cell ); 
      } 

      if( total_dangling_values == count_dangling_values() )
        try_further_possible_values_eliminaton( examined_cell );

      if( total_dangling_values == count_dangling_values() ) 
        examine_for_pairs( examined_cell );
      
      if( total_dangling_values == count_dangling_values() ) 
        examine_for_triplets( examined_cell );
 
      if(examined_cell-> return_cell_value() == 0) // the final value of examined sudoku cell is still unknown  
        examined_unsolved_cells_queue.push( examined_cell );
      
      unsolved_cells_queue.pop(); // emptying the prioroty queue 
    }  
    // Moving all the examined (and still, without a final value) cells back to unsolved_cells_queue.  
    unsolved_cells_queue = std::move( examined_unsolved_cells_queue ); 
    cout << "still "<< unsolved_cells_queue.size() << " cells to do"<<endl;

    print_grid();
  }

  // Last, sanity cheching on the result
  [&]() -> void {  
   
    set<int> collected_values; 
    for(size_t i=0; i<9; i++) {
      for(size_t j=0; j<9; j++) collected_values.insert(sudoku_array[i][j]-> return_cell_value());
      
      if(collected_values != all_possible_values) cout <<"Problem on row "<<i<<", bad result"<<endl; 
      collected_values.clear(); 
    }
   
    collected_values.clear(); 
    for(size_t i=0; i<9; i++) {
      for(size_t j=0; j<9; j++) collected_values.insert(sudoku_array[j][i]-> return_cell_value());
      
      if(collected_values != all_possible_values) cout <<"Problem on column "<<i<<", bad result"<<endl; 
      collected_values.clear(); 
    }

    collected_values.clear(); 
    for(size_t i=0; i<9; i=i+3) 
      for(size_t j=0; j<9; j=j+3) {

        size_t subblock_row = i, subblock_column = j;
        for(size_t k = subblock_row; k<subblock_row+3; k++) 
          for(size_t l = subblock_column; l<subblock_column+3; l++) 
            collected_values.insert(sudoku_array[k][l]-> return_cell_value());
      
        if(collected_values != all_possible_values) cout <<"Problem on subblock "<<i<<","<<j<<", bad result"<<endl; 
        collected_values.clear();
      }        
  }();  
  
}

// pass array cell_inputs by reference
void sudoku_solver::collect_input( char* input_file ) { 

  ifstream my_filestream;
  string fetched_line = "";
  size_t input_count = 0, i = 0;

  my_filestream.open( input_file, ios::in );
  while( std::getline( my_filestream, fetched_line ) ) {

    //Each input line must be 9 chars exactly
    if( fetched_line.size() != 9 ) { cout << "Malformed input"<<endl; exit(-1); } 
    for( size_t j=0; j<fetched_line.size(); j++ ) {

      if( isdigit(fetched_line[j]) ) cell_inputs[i][j] = fetched_line[j] - '0'; // ugly hack! 
      else cell_inputs[i][j] = -1;     
    }
    input_count += 9;
    i++;
  }
  // Missing characters
  if( input_count != 9*9 ) { cout << "Malformed input"<<endl; exit(-1); } 

  my_filestream.close();
  cout << endl;

  // If everything is ok with the input file, allocate memory for the sudoku cells 
  for( size_t i=0; i<9; i++ )
    for(size_t j=0; j<9; j++ ) {
       sudoku_array[i][j] = new sudoku_cell( cell_inputs[i][j], i, j ); 
       // cells with unknown value are placed in the unsolved_cells_queue priority queue. 
       if( cell_inputs[i][j] == -1 ) 
         unsolved_cells_queue.push( sudoku_array[i][j] ); 
    }   
}


int main( int argc, char* argv[] ) {

  sudoku_solver s_sol;
  s_sol.collect_input( argv[1] );
  s_sol.solve();
  
  return 0;
}


