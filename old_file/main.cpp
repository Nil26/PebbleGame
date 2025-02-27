//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// INCLUSIONS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//#include "stdafx.h"				// Visual studio seems to require this
#include <iostream>				// Input and output directly to screen
#include <stack>				// The stack data structure
#include <stdlib.h>				// Standard functions, including random number generation
#include <math.h>				// Basic math functions
#include <vector>				// C++ vectors, variable-length lists that can (unlike stacks and queues) have their middle elements accessed
#include <time.h>				// Get clock time, to measure total run time
#include <fstream>				// Input and output to files
//#include <string>
//#include <sstream>
#include <iomanip>

using namespace std;			// Use namespace so you can type cout instead of std::cout



 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// DECLARATIONS OF VARIABLES

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






const int ll = 64;																					// The number of vertices on a side of the lattice
const float pcrit = (float)(2/3);											// The probability of cross bracing so that the mean number of vertices is isostatic
//const int ncrit = 2 * ll - 3;																		// The isostatic number of crossbraces
const int size = ll*ll;																				// The number of vertices in the graph
short pc[size];				// Creates the pebble count at each vertex.
short occ[size];             // Says whether the site is occupied with a particle

short rcluster[size];       // Says whether the site is in a rigid cluster

const int nmax = 3 * ll;    // apparently used to keep track of number of floppy modes
int nlist[nmax];            // apparently used to keep track of number of floppy modes?
int numparts = 0;           // the number of particles (not pebbles) present in the system
int trialcount = 0;             // which trial you are on in a multitrial run
int numtries;               // within a trial, how many times you've tried to add a cluster.
float globalc;              // the clustering coefficient, set to global to increase scope
FILE *file;

#define EMPTY -size-1

ofstream myfile;																					// The file stream




vector<int> thegraph[size];		// thegraph is the graph of all loaded edges
vector<int> rgraph[size];		// rgraph is the graph of redundant bonds that don't take up any edges
vector<int> blankgraph[size];	// blankgraph is just an array of empty vectors to reset the other graphs (probably not the best way)

//vector<int> thegraphcopy[size];     //thegraphcopy is the graph of all loaded edges to implement the DFS for finding the rigid cluster

vector<int> giantrigidcluster[size];    //giantrigidcluster is the graph for the giant rigid cluster

int numbonds;					// The number of non-redundant bonds (original bonds and crossbraces) in the system
int rbonds;						// The number of redundant bonds in the system

int giantsize;                  // The size of the giant rigid cluster
int giantindex;                 // The index for the giant rigid cluster

stack<int> placesbeen;			// The list of places been while looking for a pebble

int rigidcluster();
int giantrclustersize(int giantindex_in);
bool spanningrcluster(int giantindex_in);





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ADDING, REMOVING, AND CHECKING FOR EDGES

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// addedge adds an edge pointing from vertex i to vertex j in thegraph
void addedge(int i, int j)
{
	thegraph[i].push_back(j);
	//numbonds++;
}

// addredundant adds an edge between i and j in rgraph, the separate graph of redundant edges
void addredundant(int i, int j)
{
	rgraph[i].push_back(j);
	//rbonds++;
}


// badremoveedge removes an edge pointing from i to j in the graph.
// I call it "bad" because it looks through all of i's elements instead of
// using which one might be open in a tree search.
// This should return an error if there is no edge pointing from i to j.
int badremoveedge(int i, int j)
{
	for (int k = 0; k < thegraph[i].size(); k++)
	{
		if (thegraph[i].at(k) == j)
		{
			thegraph[i].erase(thegraph[i].begin() + k);
			return 0;
		}
	}

	cout << "I tried to remove an edge pointing from " << i << " to " << j << " but I couldn't find one.\n";

	return 0;
}


// contains returns 1 if there is a non-redundant edge pointing from i to j (but doesn't check j to i) and 0 otherwise
bool contains(int i, int j)
{
	for (int k = 0; k < thegraph[i].size(); k++)
	{
		if (thegraph[i].at(k) == j)
		{
			return 1;
		}
	}


	return 0;
}


// rcontains returns 1 if there is a redundant edge pointing from i to j (but doesn't check j to i) and 0 otherwise
bool rcontains(int i, int j)
{
	for (int k = 0; k < rgraph[i].size(); k++)
	{
		if (rgraph[i].at(k) == j)
		{
			return 1;
		}
	}


	return 0;
}


// isempty returns 0 if there is any kind of redundant or nonredundant brace pointing in either direction between i and j, and 1 otherwise
// Modification required: isempty only looks at one type of bond, not the six (or at least, three) kinds of the triangular lattice
bool isempty(int i, int j)
{
	if(contains(i, j) || contains(j, i) || rcontains(i, j) || rcontains(j, i))
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

bool isempty(int i)
{
	return isempty(i, i + ll + 1);
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// FINDING NEIGHBORS ON THE TRIANGULAR LATTICE WITH PBC

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int moveright(int site)
{
    if(site % ll == ll - 1)
    {
        return site - ll + 1;
    }
    else
    {
        return site + 1;
    }
}

int moveleft(int site)
{
    if(site % ll == 0)
    {
        return site + ll - 1;
    }
    else
    {
        return site - 1;
    }
}

int moveup(int site)
{
    if(site >= ll*ll-ll)
    {
        return site + ll - ll * ll;
    }
    else
    {
        return site + ll;
    }
}

int movedown(int site)
{
    if(site <= ll - 1)
    {
        return site - ll + ll * ll;
    }
    else
    {
        return site - ll;
    }
}


// moving to the right
int dir1(int site)
{
    return moveright(site);
}


// moving up (geometrically, up and to the right)
int dir2(int site)
{
    return moveup(site);
}


// moving up and to the left (the moves should commute)
int dir3(int site)
{
    return moveleft(moveup(site));
}

int dir4(int site)
{
    return moveleft(site);
}

int dir5(int site)
{
    return movedown(site);
}

int dir6(int site)
{
    return moveright(movedown(site));
}

int choosedir(int site, int d)
{
    switch(d)
    {
    case 1:
        return dir1(site);
    case 2:
        return dir2(site);
    case 3:
        return dir3(site);
    case 4:
        return dir4(site);
    case 5:
        return dir5(site);
    case 6:
        return dir6(site);
    }
    
    return false;

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// RANDOM FUNCTIONS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// randprob() just returns a random number between zero and one.
float randprob()
{
    float f = rand();
    return f / RAND_MAX;
}



// randsite0 finds a random plaquette that may or may not be occupied
int randsite0()
{
    return (rand() % size);
}

// randsite0(skip) finds a random plaquette that skips the first and last skip rows and columns
int randsite0(int skip)
{
	return (skip+ rand() % (ll - 2*skip) + ll * (skip+(rand() % (ll - 2*skip))));
}

// randsite() finds a random plaquette that is not occupied
int randsite()
{
	int sitetoreturn = 0;
	bool empty = 0;
	int bailoutcounter = 0;
	while (empty == 0 && bailoutcounter < 10 * ll*ll)
	{
		bailoutcounter++;
		sitetoreturn = randsite0();
		empty = isempty(sitetoreturn);
	}

	if (bailoutcounter == 10 * ll *ll)
	{
		cout << "ERROR: I couldn't find an empty site! Bailout counter is: " << bailoutcounter << endl;
	}

	return sitetoreturn;
}

// randsite(skip) finds a random unoccupied site that is not in the first or last skip rows or columns
int randsite(int skip)
{
	int sitetoreturn = 0;
	bool empty = 0;
	int bailoutcounter = 0;
	while (empty == 0 && bailoutcounter < 10 * ll*ll)
	{
		bailoutcounter++;
		sitetoreturn = randsite0(skip);
		empty = isempty(sitetoreturn);
	}

	if (bailoutcounter == 10 * ll *ll)
	{
		cout << "ERROR: I couldn't find an empty site! Bailout counter is: " << bailoutcounter << endl;
	}

	return sitetoreturn;
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ADDING BONDS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




// This takes places been and reverses all the edges on the path.
// So if placesbeen is from 2 to 12 to 5
// It should remove edges from 2 to 12 and from 12 to 5
// and add edges from 5 to 12 and 12 to 2
void reversepath()
{
	int starter = placesbeen.top();				// We start at the last place in the path, the site where we found a pebble
	pc[starter]--;								// We remove a pebble from this site

	placesbeen.pop();							// We remove this site from our path, but it is still stored in starter

	int ender;



	while (placesbeen.size()>0)					// While there are still sites in the path...
	{
		ender = placesbeen.top();				// ender becomes the most recent site
		addedge(starter, ender);				// and we point an edge from starter to ender
		badremoveedge(ender, starter);			// and remove an edge from ender to starter
		starter = ender;						// then ender becomes the new starter
		placesbeen.pop();						// and is removed from the path


	}


	pc[ender]++; // Finally, we add a pebble to the first site in the path, we've moved a pebble from the end to the start, reversing edges along the way
}


// findpebble executes a depth first search for pebbles in the lattice, starting at i
// and searching all the other vertices
// placesbeen must be empty before findpebble is called
// if findpebble finds a pebble, it should set placesbeen to a path from i to the site with the pebble
// Otherwise, it should leave placesbeen blank.
// It returns 1 if a pebble was found and 0 if it wasn't or if the path wasn't empty to start
bool findpebble(int i)
{
	if (placesbeen.size() > 0)						// We shouldn't ever call findpebble unless our path has been cleared, or something didn't close like it should
	{
		cout << "I tried to find a pebble, but when I started the path stack wasn't empty, it had size " << placesbeen.size() << endl;
		return 0;
	}
	else
	{
		bool beenthere[size] = {};			// We create an array to say whether we've been to these sites before (or should we create a global array and just set it to zero here?)
		beenthere[i] = 1;					// which we haven't, to start, except for our starting site (and any skip sites)
		placesbeen.push(i);					// start our path at the starting site
		// cout << "Current location: " << placesbeen.top() << endl;

		int cl;								// cl is the current location
		int prosp;							// prosp is a vertex connected to cl that we might prospectively move to
		//cout << "Current location is " << placesbeen.top() << endl;

		while (placesbeen.size() > 0)										// Until we are forced to retreat all the way back to the first vertex...
		{
			cl = placesbeen.top();											// Our current location is the last place in the path
			// cout << "From the top, our current location is " << cl << endl;
			for (int index1 = 0; index1 < thegraph[cl].size(); index1++) // for each place we can go from our current location...
			{
				prosp = thegraph[cl].at(index1);							// Our prospective location is one of the places we can go to from cl
				//cout << "The prospective vertex we consider is " << prosp << endl;
				if (beenthere[prosp] == 0)									// if we haven't been there before...
				{
					// cl = prosp; // move our current location to there
					//cout << "We are moving to " << prosp << endl;

					placesbeen.push(prosp);									// Then we add it to the path
					//cout << "Current location after pushing: " << placesbeen.top() << endl;


					if (pc[placesbeen.top()] > 0){ return 1; }				// If our new site has a pebble, quit looking for pebbles and say we found one

					beenthere[prosp] = 1;									// Otherwise mark it as having been visited, but keep looking for a pebble
					//cout << "Current location is, in the for loop " << placesbeen.top() << endl;

					break;													// And break the for loop-- no point in continuing to explore cl's neighbors
																			// Once we've found one to move onto

				}

			}


			if (cl == placesbeen.top() || placesbeen.size() == 0) // If, after the for loop, the new top of the path is the same as the old one,
			{							// Then we didn't move anywhere, so we need to pop off the last vertex and retreat
				placesbeen.pop();

				/*
				if (placesbeen.size() == 0){ cout << "The stack is empty." << endl; }
				else
				{
				cout << "Current location after popping: " << placesbeen.top() << endl;
				}
				*/
			}



		}

		return 0;


	}
}


// When called with two arguments, findpebble skips over the second site to avoid infinitely swapping pebble
// between the two sites the brace connects, by marking skip as a place that we've already been
bool findpebble(int i, int skip)
{
	if (placesbeen.size() > 0)
	{
		cout << "I tried to find a pebble, but when I started the path stack wasn't empty, it had size " << placesbeen.size() << endl;
		return 0;
	}
	else
	{
		bool beenthere[size] = {};			// We create an array to say whether we've been to these sites before
		beenthere[i] = 1;					// which we haven't, to start, except for our starting site (and any skip sites)
		beenthere[skip] = 1;
		placesbeen.push(i);					// start our path at the starting site
		// cout << "Current location: " << placesbeen.top() << endl;


		int cl;								// cl is the current location
		int prosp;							// prosp is a vertex connected to cl that we might prospectively move to
		//cout << "Current location is " << placesbeen.top() << endl;

		while (placesbeen.size() > 0) // Until we are forced to retreat all the way back to the first vertex...
		{
			cl = placesbeen.top();
			// cout << "From the top, our current location is " << cl << endl;
			for (int index1 = 0; index1 < thegraph[cl].size(); index1++) // for each place we can go from our current location...
			{
				prosp = thegraph[cl].at(index1);
				//cout << "The prospective vertex we consider is " << prosp << endl;
				if (beenthere[prosp] == 0) // if we haven't been there before...
				{
					// cl = prosp; // move our current location to there
					//cout << "We are moving to " << prosp << endl;

					placesbeen.push(prosp);	  // and add it to the path
					//cout << "Current location after pushing: " << placesbeen.top() << endl;
					if (pc[placesbeen.top()] > 0){ return 1; }		// If our new site has a pebble, quit looking for pebbles
					beenthere[prosp] = 1;
					//cout << "Current location is, in the for loop " << placesbeen.top() << endl;
					break;


					/*
					if (pc[cl]>0)			// Then, if there is a pebble at our new location, stop looking for pebbles
					{
					return 0;
					}
					else
					{
					break;					  // Otherwise, just quit moving through the old location's neighbors
					}
					*/
				}

			}


			if (cl == placesbeen.top() || placesbeen.size() == 0) // If, after the for loop, the new top of the path is the same as the old one,
			{							// Then we didn't move anywhere, so we need to pop off the last vertex and retreat
				placesbeen.pop();

				/*
				if (placesbeen.size() == 0){ cout << "The stack is empty." << endl; }
				else
				{
				cout << "Current location after popping: " << placesbeen.top() << endl;
				}
				*/
			}



		}

		return 0;


	}
}


// loadsite looks for pebbles and moves them onto i until i has two pebbles or it stops finding pebbles
bool loadsite(int i)
{
	if (placesbeen.size() == 0)
	{
		while (pc[i] < 2 && findpebble(i)) // while the site is not loaded and you are finding pebbles...
			// c++ documentation says && short circuits, so you shouldn't
			// even look for a pebble if the site is loaded.
		{
			reversepath();					// reverse the path, which also shifts the pebbles
		}
	}
	else
	{
		cout << "I tried to find a pebble, but when I started the path stack wasn't empty, it had size " << placesbeen.size() << endl;
	}

	if (pc[i] == 2) { return 1; }		// return 1 if the site loaded successfully.
	else{ return 0; }
}

// loadsite with two arguments does the same thing, but won't try to take pebbles from skip to move onto i
// Which is to keep the two sites from swapping the three pebbles between themselves endlessly
bool loadsite(int i, int skip)
{
	if (placesbeen.size() == 0)
	{
		while (pc[i] < 2 && findpebble(i, skip)) // while the site is not loaded and you are finding pebbles...
			// C++ documentation says && short-circuits, so this shouldn't even
			// look for a pebble if the site is loaded.
		{
			reversepath();					// reverse the path, which also shifts the pebbles
		}
	}
	else
	{
		cout << "I tried to load site " << i << " but the placesbeen stack wasn't empty.\n";
	}
	if (pc[i] == 2) { return 1; }		// return 1 if the site loaded successfully.
	else{ return 0; }
}

// loadsites tries to move pebbles until there are two on both sites i and j
bool loadsites(int i, int j)
{

	while (pc[j] < 2 && findpebble(j))
	{
		reversepath();
	}

	while (pc[i] < 2 && findpebble(i, j))	// while the first site is not loaded and you are finding pebbles skipping
	{										// second site, load the first site, then load the second.
		reversepath();
		while (pc[j] < 2 && findpebble(j))
		{
			reversepath();
		}
	}
	if (pc[i] == 2 && pc[j] == 2) { return 1; }
	else{ return 0; }
}

// addbond tries to load the sites. If it succeeds, it adds an edge from i to j and takes a pebble from i. Otherwise, it adds a redundant edge
void addbond(int i, int j)
{
	if (numbonds < 2 * size - 3 && loadsites(i, j))			// If there are at least four pebbles left, we try to load the sites
	{
		addedge(i, j);				// if we succeed, we add the edge from i to j and remove a pebble from i
		numbonds++;
		pc[i]--;
	}
	else
	{
		addredundant(i, j);			// otherwise, we leave the pebbles where we shuffled them and
		rbonds++;
	}								// place only a redundant bond



}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// REPORTING TO SCREEN

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void printfilename(float cval, int tval)
{
    char mychar[] = "cxxxtxxxx.txt";
    int start =sizeof(mychar)-13;
    cval = cval * 100;
    int intc = cval;
    mychar[start+4]='0'+ (tval % 10000)/1000;
    mychar[start+5]='0'+ (tval % 1000)/100;
    mychar[start+6]='0'+ (tval % 100)/10;
    mychar[start+7]='0'+ (tval % 10);

    mychar[start]='0'+ (intc % 1000)/100;
    mychar[start+1]='0'+ (intc % 100)/10;
    mychar[start+2]='0'+ (intc % 10);


    cout << "filename: " << mychar << endl;

}

void setfilestream(float cval, int tval)
{


    char mychar[] = "data1/cxxxtxxxx.txt";
    int start =sizeof(mychar)-13;
    cval = cval * 100;
    int intc = cval;
    mychar[start+4]='0'+ (tval % 10000)/1000;
    mychar[start+5]='0'+ (tval % 1000)/100;
    mychar[start+6]='0'+ (tval % 100)/10;
    mychar[start+7]='0'+ (tval % 10);

    mychar[start]='0'+ (intc % 1000)/100;
    mychar[start+1]='0'+ (intc % 100)/10;
    mychar[start+2]='0'+ (intc % 10);

    myfile.close();
    myfile.open(mychar);

    //cout << mychar;
}



void log(int span=-1)
{
    myfile << trialcount << "\t" << ll << "\t" << globalc << "\t" << numparts << "\t" << numbonds << "\t" << rbonds << "\t" << giantindex << "\t" << giantsize <<"\t"<<span<< "\n";
}

void log2()
{
    char ss[] = "hello world";
    fputs(ss,file);
}

/*
std::string mof(int index) {
   std::ostringstream ss;
   ss << "testfile_" << index << ".txt";
   return ss.str();
}
*/


/*
char tochar(string s)
{
    char a[2];
    a[0]='h';
    a[1]='w';

    return a;
}
*/

void plot()
{
    cout << "\n\n";

    for(int j = 0; j< ll; j++)
    {
        if (j%2==0) {
            cout<<" ";
        }
        for(int i=0;i<ll;i++)
            {
                cout << setw(1) << occ[i+j*ll] << " ";
            }
    cout << endl;

    }

}

// listedges lists the edges from site i
void listedges(int i)
{
	cout << "\nvertex number " << i << " points towards the following vertices:\n";

	for (int index = 0; index < thegraph[i].size(); index++)
	{
		cout << thegraph[i].at(index) << " ";
	}

	cout << endl;
}


// listalledges lists all the edges from the sites
void listalledges()
{
	for (int index = 0; index < size; index++)
	{
		listedges(index);
	}
}


// report lists all the edges from each site and all the pebbles on each site, as well as the total number of bonds and pebbles
void report()
{
	int psum = 0;
	for (int i = 0; i < size; i++)
	{
		psum += pc[i];
	}

	cout << "The number of (nonredundant) bonds is " << numbonds << endl;
	cout << "The number or redundant bonds is " << rbonds << endl;
	cout << "The number of remaining pebbles is " << psum << "\n\n";

	for (int index = 0; index < size; index++)
	{
		cout << "vertex number " << index << ":\n";
		cout << "pebble count: " << pc[index] << endl;
		cout << "Neighbors: ";

		for (int index2 = 0; index2 < thegraph[index].size(); index2++)
		{
			cout << thegraph[index].at(index2) << " ";
		}

		cout << "\n\n";
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// FIXED-PROBABILITY TRIALS AND ATTENDANT FUNCTIONS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// initializegraph() updates numbonds, rbonds, thegraph, rgraph, placesbeen to a square graph with no braces
void initializegraph()
{
	numbonds = 2 * ll*ll - 2 * ll;
	rbonds = 0;

	for (int pcindex = 0; pcindex < size; pcindex++) // Just setting the pebble count to 2 everywhere.
	{
		pc[pcindex] = 0;
	}

	/*pc[ll*ll-1] = 2;										// Then put the 2 L pebbles on the edges

	for (int pcindex = 1; pcindex < ll; pcindex++)
	{
	pc[pcindex] = 1;
	pc[ll*pcindex] = 1;
	}
	*/

	while (placesbeen.size()>0)					// Clear the places been stack
	{
		placesbeen.pop();
	}

	for (int bondindex = 0; bondindex < size; bondindex++)
	{
		rgraph[bondindex].clear();
		thegraph[bondindex].clear();

		if (bondindex % ll != ll - 1)
		{
			thegraph[bondindex].push_back(bondindex + 1);
		}
		else
		{
			pc[bondindex]++;
		}

		if (bondindex < ll*(ll - 1))
		{
			thegraph[bondindex].push_back(bondindex + ll);
		}
		else
		{
			pc[bondindex]++;
		}
	}


}

// inittrigraph() updates numbonds, rbonds, thegraph, rgraph, placesbeen to a triangular graph with particles but no bonds
int inittrigraph(float p)
{
	numbonds = 0; // the number of bonds
    numparts = 0; // the number of particles
	rbonds = 0;   // the number of redundant bonds



	for (int pcindex = 0; pcindex < size; pcindex++) // Just setting the pebble count to 2 everywhere.
	{                                               // and setting which sites are occupied
		pc[pcindex] = 2;



		if(randprob()<p)
        {
            occ[pcindex]=1;
            numparts++;
        }
        else
        {
            occ[pcindex]=0;
        }

	}



	while (placesbeen.size()>0)					// Clear the places been stack
	{
		placesbeen.pop();
	}

	for (int bondindex = 0; bondindex < size; bondindex++)  // Clear the graphs of redundant and nonredundant bonds
	{
		rgraph[bondindex].clear();
		thegraph[bondindex].clear();

	}


    return numparts;

}



// inittrigraph() updates numbonds, rbonds, thegraph, rgraph, placesbeen to a triangular graph with particles but no bonds
int inittrigraph(float p, float c)
{
	numbonds = 0; // the number of bonds
    numparts = 0; // the number of particles
	rbonds = 0;   // the number of redundant bonds

    float p0 = p / (1 + 6 * c);



	for (int pcindex = 0; pcindex < size; pcindex++) // Just setting the pebble count to 2 everywhere.
	{                                               // and setting which sites are occupied
		pc[pcindex] = 2;
		occ[pcindex]=0;
	}



	for (int pcindex = 0; pcindex < size; pcindex++) // Just setting the pebble count to 2 everywhere.
	{                                               // and setting which sites are occupied



		if(randprob()<p0)
        {
            if(occ[pcindex]==0)
            {
                occ[pcindex]=1;
                numparts++;
            }



                if( randprob()<c && occ[dir1(pcindex)]==0)
                {
                    occ[dir1(pcindex)]=1;
                    numparts++;
                }

                if( randprob()<c && occ[dir2(pcindex)]==0)
                {
                    occ[dir2(pcindex)]=1;
                    numparts++;
                }

                if( randprob()<c && occ[dir3(pcindex)]==0)
                {
                    occ[dir3(pcindex)]=1;
                    numparts++;
                }


                if( randprob()<c && occ[dir4(pcindex)]==0)
                {
                    occ[dir4(pcindex)]=1;
                    numparts++;
                }



                if( randprob()<c && occ[dir5(pcindex)]==0)
                {
                    occ[dir5(pcindex)]=1;
                    numparts++;
                }



                if( randprob()<c && occ[dir6(pcindex)]==0)
                {
                    occ[dir6(pcindex)]=1;
                    numparts++;
                }




        }


	}



	while (placesbeen.size()>0)					// Clear the places been stack
	{
		placesbeen.pop();
	}

	for (int bondindex = 0; bondindex < size; bondindex++)  // Clear the graphs of redundant and nonredundant bonds
	{
		rgraph[bondindex].clear();
		thegraph[bondindex].clear();

	}


    return numparts;

}

// initgiantrigidcluster() initializes the giantrigidcluster graph
void initgiantrigidcluster()
{
    for (int bondindex = 0; bondindex < size; bondindex++)  // Clear the graphs of giantrigidcluster
    {
        giantrigidcluster[bondindex].clear();
        
    }
}


// initemptytrigraph() updates numbonds, rbonds, thegraph, rgraph, placesbeen to a triangular graph with no particles or bonds
int initemptytrigraph()
{
	numbonds = 0; // the number of bonds
    numparts = 0; // the number of particles
	rbonds = 0;   // the number of redundant bonds




	for (int pcindex = 0; pcindex < size; pcindex++) // Just setting the pebble count to 2 everywhere.
	{                                               // and setting which sites are occupied
		pc[pcindex] = 2;
		occ[pcindex]=0;
	}


	while (placesbeen.size()>0)					// Clear the places been stack
	{
		placesbeen.pop();
	}

	for (int bondindex = 0; bondindex < size; bondindex++)  // Clear the graphs of redundant and nonredundant bonds
	{
		rgraph[bondindex].clear();
		thegraph[bondindex].clear();
        giantrigidcluster[bondindex].clear();

	}
    
    initgiantrigidcluster();


    return numparts;

}


// addrandombraces adds a random brace to each of the (ll-1)^2 sites with probability p
void addrandombraces(float p)
{
	for (int siteindex = 0; siteindex < ll*(ll - 1); siteindex++)
	{
		if (siteindex % ll != ll - 1 && randprob() < p)
		{
			addbond(siteindex, siteindex + ll + 1);
		}
	}
}

// addtribraces adds a random brace to each of the sites with probability p, only if both are occupied
void addtribraces()
{
    int debugcount = 0;
    int newsite;
	for (int siteindex = 0; siteindex < ll*ll; siteindex++)
	{
	    if(siteindex % ll == ll-1)
        {
            newsite = siteindex + 1 - ll;
        }
        else
        {
            newsite = siteindex + 1;
        }
		if (occ[siteindex] == 1 && occ[newsite]==1)
		{
			addbond(siteindex, newsite);
			debugcount++;
		}

        if(siteindex >=ll*(ll-1))
        {
            newsite = siteindex + ll - ll * ll;
        }
        else
        {
            newsite = siteindex + ll;
        }
		if (occ[siteindex] == 1 && occ[newsite]==1)
		{
			addbond(siteindex, newsite);
			debugcount++;
		}

        if(siteindex >=ll*(ll-1))
        {
            newsite = siteindex + ll - ll * ll;
        }
        else
        {
            newsite = siteindex + ll;
        }
        if(newsite % ll ==0)
        {
            newsite = newsite + ll -1;
        }
        else
        {
            newsite = newsite - 1;
        }
		if (occ[siteindex] == 1 && occ[newsite]==1)
		{
			addbond(siteindex, newsite);
			debugcount++;
		}
	}
	//cout << "Number of added braces: " << debugcount << endl;
}


void addtripart(int site)
{
    int newsite;

    if(occ[site]==1)
    {
        // nothing
    }
    else
    {

        occ[site]=1;

        numparts++;
        for(int k = 1; k<=6; k++)
        {
            newsite = choosedir(site,k);
            if(occ[newsite]==1&&isempty(site,newsite))
            {
                addbond(site, newsite);
            }
        }
    }

}

void addtricluster(int site, float c)
{
    int oldnumparts = numparts;
    int newedgesite;

    addtripart(site);

    for(int k = 1; k<=6; k++)
    {
            if(randprob()<c)
            {
                newedgesite = choosedir(site,k);

                //cout << "old site: " << site << " new site: " << newedgesite << endl;
                addtripart(newedgesite);

                //plot();
            }

    }


    if(numparts>oldnumparts)
    {
        log();
    }

}


void addtricluster2(int site, float c) // Has already added the rigidcluster function, as well as the spanning cluster
{

    if(occ[site]==0)
    {
        int newsite;
        int numneighbors = occ[choosedir(site,1)]+occ[choosedir(site,2)]+occ[choosedir(site,3)]+occ[choosedir(site,4)]+occ[choosedir(site,5)]+occ[choosedir(site,6)];
        //cout << "base is " << 1.-c << " and power is " << 6-numneighbors << ".  ";
        //cout << "probability of placement is " << pow(1.-c,6-numneighbors) << endl;

        if(randprob()<pow(1.-c,6-numneighbors))
        {
            occ[site]=1;

            numparts++;
            for(int k = 1; k<=6; k++)
            {
                newsite = choosedir(site,k);
                if(occ[newsite]==1&&isempty(site,newsite))
                {
                    addbond(site, newsite);
                }


            }
			// choose some densities for the rigid cluster
            giantindex=rigidcluster();
            giantsize=giantrclustersize(giantindex);
            
            int span=spanningrcluster(giantindex);
            log(span);

        }
    }

}

void onetritrial2(int maxout,float c)
{
    int numattempts = 0;
    initemptytrigraph();
    trialcount = 1;
    
    while(numattempts < maxout && numparts < ll*ll)
    {
        numattempts++;
        addtricluster2(randsite0(),c);
        //plot();
    }
    
    myfile.close();
}

void tritrials(int maxout,float c,int numtrials)
{


    for(trialcount = 1; trialcount<=numtrials; trialcount++)
    {
        int numattempts = 0;
        initemptytrigraph();

        while(numattempts < maxout && numparts < ll*ll)
            {
                numattempts++;
                addtricluster(randsite0(),c);
            }
    }
}

void multictrial(int maxout,float c1,float c2, float dc,int numtrials)
{
    for(globalc = c1; globalc <= c2; globalc += dc)
    {
            for(int mtc = 14; mtc<=numtrials; mtc++)
            {
                trialcount++;
                int numattempts = 0;
                initemptytrigraph();
                /*
                if(mtc = 1 && globalc > .9)
                {
                    plot();
                    cout << "numattempts: " << numattempts << "c: " << globalc << "numtrials: " << numtrials << endl;
                }
                */
                setfilestream(globalc,mtc);



                while(numattempts < maxout && numparts < ll*ll)
                    {
                        numattempts++;

                        addtricluster2(randsite0(),globalc);
                    }
            }

    }


}




// onetrial adds random braces with probability p and reports the number of floppy modes
void onetrial(float p)
{
	initializegraph();

	addrandombraces(p);

	cout << "\nWhen the number of vertices on a side is " << ll << " and the probability of a random brace being placed is \n"
		<< p << " the number of floppy modes is\n" << 2 * ll*ll - numbonds << endl;

}

void onetritrial(float p)
{
	inittrigraph(p);

	addtribraces();

	cout << "\nWhen the number of vertices on a side is " << ll << " and the probability of a random brace being placed is \n"
		<< p << " the number of floppy modes is\n" << 2 * ll*ll - numbonds << endl;

}

// multitrial adds random braces trial times and returns the mean number of floppy modes
float multitrial(float p, int trial)
{

	int numfloppy = 0;

	for (int t = 0; t < trial; t++)
	{
		initializegraph();

		addrandombraces(p);
		numfloppy +=(2 * ll *ll - numbonds);
	}

	return ((float)numfloppy / trial);

}


// multitrial adds random braces trial times and returns the mean number of floppy modes
float multitritrial(float p, int trial)
{

	float fracfloppy = 0;

	for (int t = 0; t < trial; t++)
	{
		inittrigraph(p);

		addtribraces();

		myfile << ll << "\t" << p << "\t" << numparts << "\t" << numbonds << "\t" << rbonds << "\n";

	}



	return ((float)fracfloppy / trial);

}


// multitrial adds random braces trial times and returns the mean number of floppy modes
float multitritrial(float p, float c, int trial)
{

	float fracfloppy = 0;

	for (int t = 0; t < trial; t++)
	{
		inittrigraph(p,c);

		addtribraces();

		myfile << ll << "\t" << p << "\t" << c << "\t" << numparts << "\t" << numbonds << "\t" << rbonds << "\n";

	}



	return ((float)fracfloppy / trial);

}


// prun goes from probability p1 to p2 by dp, calculating the mean number of floppy modes based on trials at each point, and writing the results to a file
void pruntri(float p1, float p2, float dp, int trials)
{

	for (float p = p1; p <= p2; p += dp)
	{


		for (int t = 0; t < trials; t++)
		{
			multitritrial(p,trials);
		}



	}
}

// prun goes from probability p1 to p2 by dp, calculating the mean number of floppy modes based on trials at each point, and writing the results to a file
void pruntri(float p1, float p2, float dp, int trials,float c)
{

	for (float p = p1; p <= p2; p += dp)
	{


		for (int t = 0; t < trials; t++)
		{
			multitritrial(p,c,trials);
		}



	}
}

// prun goes from probability p1 to p2 by dp, calculating the mean number of floppy modes based on trials at each point, and writing the results to a file
void pruntri(float p1, float p2, float dp, float c1, float c2, float dc, int trials)
{

	for (float p = p1; p <= p2; p += dp)
	{

	    for(float c = c1; c <= c2; c += dc)
	    {
            for (int t = 0; t < trials; t++)
            {
                multitritrial(p,c,trials);
            }

	    }






	}
}








//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// fixed-number trials and attendant functions

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// addrandombracesfixedn adds ntoadd random braces to sites, taking care that braces do not overlap one another.
// it should give an error if there aren't enough sites to accommodate all the braces.
void addrandombracesfixedn(int ntoadd)
{

	int nadded = 0;
	int thesite;

	while (nadded < ntoadd)
	{
		thesite = randsite();
		addbond(thesite, thesite + ll + 1);
		nadded++;
	}


}

// onetrial adds ntoadd random braces and reports the number of floppy modes
void onetrialfixedn(int ntoadd)
{
	initializegraph();

	addrandombracesfixedn(ntoadd);

	cout << "\nwhen the number of vertices on a side is " << ll << " and the number of random braces placed is \n"
		<< ntoadd << " the number of floppy modes is\n" << 2 * ll*ll - numbonds << endl;

}

// multitrialfixedn adds ntoadd random braces trial times and returns the mean number of floppy modes
float multitrialfixedn(int ntoadd, int trial)
{

	int numfloppy = 0;

	for (int t = 0; t < trial; t++)
	{
		initializegraph();

		addrandombracesfixedn(ntoadd);
		numfloppy += 2 * ll *ll - numbonds;
	}

	return ((float)numfloppy / trial);

}

void nrun(int trial)
{

	int nadded;
	int thesite;




	for (int nindex = 0; nindex < nmax; nindex++)
	{
		nlist[nindex] = 0;
	}





	for (int tindex = 0; tindex < trial; tindex++)
	{
		initializegraph();
		nadded = 0;



		while (nadded < nmax)
		{
			nlist[nadded] += 2 * size - numbonds;
			thesite = randsite();
			addbond(thesite, thesite + ll + 1);
			nadded++;
		}
	}

	for (int nval = 0; nval < nmax; nval++)
	{
		myfile << nval << "\t" << ((float)nlist[nval] / trial) << "\n";
	}

}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// rigidity testing

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// a site is rigid if a bond we placed there would be redundant, when loadsite() returns 0.
// so we just try to load the site to see if it is rigid

bool isrigidsquare(int i)
{

	if (i % ll == ll - 1 || i >= ll*(ll - 1))
	{
		cout << "i tried to load a site that fell out of the system!" << endl;
	}

	return !loadsites(i, i + ll + 1);
}


bool isrigidbulk(int skip, int tries)
{
	short yescount = 0;

	for (int tindex = 0; tindex < tries; tindex++)
	{
		if (isrigidsquare(randsite()))
		{
			yescount++;
		}
	}
	if (yescount == 0)
	{
		return 0;
	}
	else
	{
		if (yescount < tries)
		{
			cout << "i got ambiguous results when testing bulk rigidity, which shouldn't happen for large systems\n";
		}
		return 1;
	}

}

float probrigidbulk(int numbraces, int trials,int skip, int tries)
{
	int numrigid = 0;

	for (int tindex = 0; tindex < trials; tindex++)
	{
		initializegraph();
		addrandombracesfixedn(numbraces);
		numrigid += isrigidbulk(skip, tries);
	}

	return ((float)numrigid / trials);

}

void probrigidbulkrun(int numstart, int numfinish, int trials, int skip, int tries)
{


	for (int numbraces = 2 * ll - 3 + numstart; numbraces < 2 * ll - 3 + numfinish; numbraces++)
	{

		int numrigid = 0;

		for (int tindex = 0; tindex < trials; tindex++)
		{
			initializegraph();
			addrandombracesfixedn(numbraces);
			numrigid += isrigidbulk(skip, tries);
		}

		myfile <<numbraces - (2*ll-3) << "\t" << ((float)numrigid / trials) << endl;

		cout << numbraces - (2 * ll - 3)  << "\t" << ((float)numrigid / trials) << endl;

	}


}

void rigiditymaptoscreen()
{
	for (int jindex = 0; jindex < ll - 1; jindex++)
	{
		for (int iindex = 0; iindex < ll - 1; iindex++)
		{
			cout << isrigidsquare(iindex + ll*jindex) << "\t";
		}

		cout << "\n";
	}
}

void rigiditymap()
{
	for (int jindex = 0; jindex < ll - 1; jindex++)
	{
		for (int iindex = 0; iindex < ll - 1; iindex++)
		{
			if (isempty(iindex + ll*jindex))
			{
				myfile << isrigidsquare(iindex + ll*jindex) << "\t";
			}
			else
			{
				myfile << 2 << "\t";
			}
		}

		myfile << "\n";
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// rigid cluster

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initrcluster() // initialize the rcluster array to store the network
{
    for(int i=0; i<=size-1; i++)
    {
        rcluster[i]=occ[i];
    }
    giantsize=0;
}

bool isredundant(int i, int j)  //see if the test bond between (i,j) is redundant(dependent)
{
    if (numbonds < 2 * size - 3 && loadsites(i, j))			// if there are at least four pebbles left, we try to load the sites
    {
        return 0;				// if we succeed, then the edge from i to j is independent
    }
    else
    {
        return 1;			// otherwise, the edge is redundant(dependent)
    }

}

int rigidcluster() // mark the rigid clusters and put into rcluster, and return the size of the giant rigid cluster.
{
    initrcluster();
    int rcnum=1; //index of the rigid cluster
    int giantsize=2; //if there exists any bond, the smallest possible giant rigid cluster size is 2
    int giantindex_process=0; //the index for the giant rigid cluster (in this function)

    for(int i=0; i<=size-1; i++)
    {
        if (rcluster[i]==1 && !thegraph[i].empty())
        {
            for (int index=thegraph[i].size()-1; index>=0; index--)
            {
                if (rcluster[thegraph[i][index]]==1)
                { //the bond exists and not in some rigid clusters

                    int rclustersize=2;//the size of this rigid cluster
                    rcnum++;
                    rcluster[i]=rcnum;
                    rcluster[thegraph[i][index]]=rcnum;
                    for (int t=0; t<=size-1; t++)
                    {
                        if (rcluster[t]==1)
                        {   //the test site which is not in some rigid clusters
                            if (isredundant(i, t) && isredundant(thegraph[i][index], t))
                            {   //see if t is in the rigid cluster of (i,thegraph[i][index])
                                rcluster[t]=rcnum;
                                rclustersize++;
                            }
                        }
                    }

                    if (rclustersize>=giantsize) //update the size of the giant rigid cluster
                    {
                        giantsize=rclustersize;
                        giantindex_process=rcnum;
                    }

                }
            }
        }
    }

    // pick out the giant rigid cluster and store it in the vector "giantrigidcluster" (!!! we need it to become a undirected adjacent list)
    initgiantrigidcluster(); //empty the vector array first
    for (int i=0; i<=size-1; i++) {
            if (rcluster[i]==giantindex_process) {
                for (int index=thegraph[i].size()-1; index>=0; index--) {
                    if (rcluster[thegraph[i][index]]==giantindex_process) { // if the bond is connected to another site in the giant rigid cluster
                        giantrigidcluster[i].push_back(thegraph[i][index]);
                        giantrigidcluster[thegraph[i][index]].push_back(i);// The giantrigidcluster graph is a undirected graph
                    }
                }
            }
    }
    
    return giantindex_process;
}

void plotrcluster()
{
    cout << "\n\n";
    
    for(int j = 0; j< ll; j++)
    {
        if (j%2==0) {
            cout<<" ";
        }
        for(int i=0;i<ll;i++)
        {
            cout << setw(1) << rcluster[i+j*ll] << " ";
        }
        cout << endl;
        
    }
    
}

int giantrclustersize(int giantindex_in) //To find the size for the giant rigid cluster based on the index of that cluster
{
    giantsize=0;    //Initialize the giantsize varialbe
    int giantsize_return=0;
    for (int i=0; i<=size-1; i++)
    {
        if (rcluster[i]==giantindex_in)
        {
            giantsize_return++;
        }
    }
    if (giantindex_in==0 || giantindex_in==1) {
        return 1;
    } else {
        return giantsize_return;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// SPANNING RIGID CLUSTER

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// We need to pick out the giant rigid cluster from the network and then determine if it is the spanning cluster.

bool spanningrcluster(int giantindex_in)
{
    // Do the DFS in the giant rigid cluster
    int beenthere[size]={};
    std::fill_n(beenthere,size,EMPTY);
    stack<int> DFS_rcluster;
    // Find a starting point
    int v_start;
    for (v_start=0; v_start<=size-1; v_start++) {
        if (rcluster[v_start]==giantindex_in) {
            break;
        }
    }
    
    //v_start is the starting point for the DFS
    beenthere[v_start]=0;
    DFS_rcluster.push(v_start);
    int cl; //cl is the current location
    int prosp; //prosp is the prospective location we're going to move to
    while (DFS_rcluster.size()>0) // make sure the DFS is in the giant rigid cluster
    {
        cl=DFS_rcluster.top();
        for (int index1=0; index1<giantrigidcluster[cl].size(); index1++)
        {
            prosp=giantrigidcluster[cl].at(index1);
            if (beenthere[prosp]==EMPTY)
            {
                DFS_rcluster.push(prosp);
                //Do sth here
                //beenthere[prosp]=1;
                
                
                //Get the x-displacement for prosp from the x-displacement for cl (x_prosp=x_cl +1 OR +0 OR -1)
                if (prosp==choosedir(cl, 2)||prosp==choosedir(cl, 5))
                {
                    beenthere[prosp]=beenthere[cl];
                } else if (prosp==choosedir(cl, 1)||prosp==choosedir(cl, 6))
                {
                    beenthere[prosp]=beenthere[cl]+1;
                } else if (prosp==choosedir(cl, 3)||prosp==choosedir(cl, 4))
                {
                    beenthere[prosp]=beenthere[cl]-1;
                }
                else
                {
                    cout<<"Whoooops, sth goes wrong."<<endl;
                }
                
                break;
            }
            else //Existing marked displacement can be compared to current place to see if the spanning cluster exists.
            {
                if (abs(beenthere[cl]-beenthere[prosp])==ll-1)
                {
                    return true;
                }
            }
        }
        if (cl==DFS_rcluster.top()||DFS_rcluster.size()==0) {
            DFS_rcluster.pop();
        }
    }
    return false;// After the DFS still no finding about the spanning cluster, then it is not spanning.
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// MAIN FUNCTION

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
	//myfile.open("tridata2.txt");										// Opens a file to write to
	//myfile << "test\n";

	clock_t t1, t2;													// Creates clock variables to track the runtime
	t1 = clock();
	//unsigned seed = t1;
	srand(time(NULL));

	//multictrial(ll*ll*100000,0.0,0.1,.2,20); //That's it!

	float cfor = 0.0; // correlation constant
	int trial = 1; // trial counting

	setfilestream(cfor,trial);
	onetritrial2(ll*ll*300000,cfor);

	t2 = clock();
	float clocktime((float)t2 - (float)t1);
	cout << "\n The total run time was " << clocktime / CLOCKS_PER_SEC << endl;
	cin.ignore();
	return 0;
}

