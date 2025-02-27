//
// Created by Zeb & Shang.
//

#include <iostream>
#include <vector>
#include <stack>
#include <fstream>
#include <math.h>				// Basic math functions
#include "siteRP.h"

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // ADDING, REMOVING, AND CHECKING FOR EDGES

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // addedge adds an edge pointing from vertex i to vertex j in thegraph
    void SiteRP::addedge(int i, int j) {
        thegraph[i].push_back(j);
    }

    // addredundant adds an edge between i and j in rgraph, the separate graph of redundant edges
    void SiteRP::addredundant(int i, int j) {
        rgraph[i].push_back(j);
    }

    // badremoveedge removes an edge pointing from i to j in the graph.
// I call it "bad" because it looks through all of i's elements instead of
// using which one might be open in a tree search.
// This should return an error if there is no edge pointing from i to j.
    int SiteRP::badremoveedge(int i, int j) {
        for (int k = 0; k < thegraph[i].size(); k++) {
            if (thegraph[i].at(k) == j) {
                thegraph[i].erase(thegraph[i].begin() + k);
                return 0;
            }
        }
        std::cout << "I tried to remove an edge pointing from " << i << " to " << j << " but I couldn't find one.\n";
        return 0;
    }

// contains returns 1 if there is a non-redundant edge pointing from i to j (but doesn't check j to i) and 0 otherwise
    bool SiteRP::contains(int i, int j) {
        for (int k = 0; k < thegraph[i].size(); k++) {
            if (thegraph[i].at(k) == j) {
                return 1;
            }
        }
        return 0;
    }


// rcontains returns 1 if there is a redundant edge pointing from i to j (but doesn't check j to i) and 0 otherwise
    bool SiteRP::rcontains(int i, int j) {
        for (int k = 0; k < rgraph[i].size(); k++) {
            if (rgraph[i].at(k) == j) {
                return 1;
            }
        }


        return 0;
    }


// isempty returns 0 if there is any kind of redundant or nonredundant brace pointing in either direction between i and j, and 1 otherwise
// Modification required: isempty only looks at one type of bond, not the six (or at least, three) kinds of the triangular lattice
    bool SiteRP::isempty(int i, int j) {
        if (contains(i, j) || contains(j, i) || rcontains(i, j) || rcontains(j, i)) {
            return 0;
        }
        else {
            return 1;
        }
    }

    bool SiteRP::isempty(int i) {
        return isempty(i, i + ll + 1);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// FINDING NEIGHBORS ON THE TRIANGULAR LATTICE WITH PBC

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    int SiteRP::moveright(int site) {
        if (site % ll == ll - 1) {
            return site - ll + 1;
        }
        else {
            return site + 1;
        }
    }

    int SiteRP::moveleft(int site) {
        if (site % ll == 0) {
            return site + ll - 1;
        }
        else {
            return site - 1;
        }
    }

    int SiteRP::moveup(int site) {
        if (site >= ll * ll - ll) {
            return site + ll - ll * ll;
        }
        else {
            return site + ll;
        }
    }

    int SiteRP::movedown(int site) {
        if (site <= ll - 1) {
            return site - ll + ll * ll;
        }
        else {
            return site - ll;
        }
    }


// moving to the right
    int SiteRP::dir1(int site) {
        return moveright(site);
    }


// moving up (geometrically, up and to the right)
    int SiteRP::dir2(int site) {
        return moveup(site);
    }


// moving up and to the left (the moves should commute)
    int SiteRP::dir3(int site) {
        return moveleft(moveup(site));
    }

    int SiteRP::dir4(int site) {
        return moveleft(site);
    }

    int SiteRP::dir5(int site) {
        return movedown(site);
    }

    int SiteRP::dir6(int site) {
        return moveright(movedown(site));
    }

    int SiteRP::choosedir(int site, int d) {
        switch (d) {
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
    float SiteRP::randprob() {
        float f = rand();
        return f / RAND_MAX;
    }

// randsite0 finds a random plaquette that may or may not be occupied
    int SiteRP::randsite0() {
        return (rand() % size);
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
    void SiteRP::reversepath() {
        int starter = placesbeen.top();                // We start at the last place in the path, the site where we found a pebble
        pc[starter]--;                                // We remove a pebble from this site

        placesbeen.pop();                            // We remove this site from our path, but it is still stored in starter

        int ender;

        while (placesbeen.size() > 0)                    // While there are still sites in the path...
        {
            ender = placesbeen.top();                // ender becomes the most recent site
            addedge(starter, ender);                // and we point an edge from starter to ender
            badremoveedge(ender, starter);            // and remove an edge from ender to starter
            starter = ender;                        // then ender becomes the new starter
            placesbeen.pop();                        // and is removed from the path
        }
        pc[ender]++; // Finally, we add a pebble to the first site in the path, we've moved a pebble from the end to the start, reversing edges along the way
    }


// findpebble executes a depth first search for pebbles in the lattice, starting at i
// and searching all the other vertices
// placesbeen must be empty before findpebble is called
// if findpebble finds a pebble, it should set placesbeen to a path from i to the site with the pebble
// Otherwise, it should leave placesbeen blank.
// It returns 1 if a pebble was found and 0 if it wasn't or if the path wasn't empty to start
    bool SiteRP::findpebble(int i) {
        if (placesbeen.size() > 0)                        // We shouldn't ever call findpebble unless our path has been cleared, or something didn't close like it should
        {
            std::cout << "I tried to find a pebble, but when I started the path stack wasn't empty, it had size " <<
            placesbeen.size() << std::endl;
            return 0;
        }
        else {
            bool beenthere[size] = {};            // We create an array to say whether we've been to these sites before (or should we create a global array and just set it to zero here?)
            beenthere[i] = 1;                    // which we haven't, to start, except for our starting site (and any skip sites)
            placesbeen.push(i);                    // start our path at the starting site
            // cout << "Current location: " << placesbeen.top() << endl;

            int cl;                                // cl is the current location
            int prosp;                            // prosp is a vertex connected to cl that we might prospectively move to
            //cout << "Current location is " << placesbeen.top() << endl;

            while (placesbeen.size() > 0)                                        // Until we are forced to retreat all the way back to the first vertex...
            {
                cl = placesbeen.top();                                            // Our current location is the last place in the path
                // cout << "From the top, our current location is " << cl << endl;
                for (int index1 = 0;
                     index1 < thegraph[cl].size(); index1++) // for each place we can go from our current location...
                {
                    prosp = thegraph[cl].at(index1);                            // Our prospective location is one of the places we can go to from cl
                    //cout << "The prospective vertex we consider is " << prosp << endl;
                    if (beenthere[prosp] == 0)                                    // if we haven't been there before...
                    {
                        // cl = prosp; // move our current location to there
                        //cout << "We are moving to " << prosp << endl;

                        placesbeen.push(prosp);                                    // Then we add it to the path
                        //cout << "Current location after pushing: " << placesbeen.top() << endl;


                        if (pc[placesbeen.top()] > 0) { return 1; }                // If our new site has a pebble, quit looking for pebbles and say we found one

                        beenthere[prosp] = 1;                                    // Otherwise mark it as having been visited, but keep looking for a pebble
                        //cout << "Current location is, in the for loop " << placesbeen.top() << endl;

                        break;                                                    // And break the for loop-- no point in continuing to explore cl's neighbors
                        // Once we've found one to move onto
                    }
                }


                if (cl == placesbeen.top() || placesbeen.size() == 0) // If, after the for loop, the new top of the path is the same as the old one,
                {                            // Then we didn't move anywhere, so we need to pop off the last vertex and retreat
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
    bool SiteRP::findpebble(int i, int skip) {
        if (placesbeen.size() > 0) {
            std::cout << "I tried to find a pebble, but when I started the path stack wasn't empty, it had size " <<
            placesbeen.size() << std::endl;
            return 0;
        }
        else {
            bool beenthere[size] = {};            // We create an array to say whether we've been to these sites before
            beenthere[i] = 1;                    // which we haven't, to start, except for our starting site (and any skip sites)
            beenthere[skip] = 1;
            placesbeen.push(i);                    // start our path at the starting site
            // cout << "Current location: " << placesbeen.top() << endl;


            int cl;                                // cl is the current location
            int prosp;                            // prosp is a vertex connected to cl that we might prospectively move to
            //cout << "Current location is " << placesbeen.top() << endl;

            while (placesbeen.size() > 0) // Until we are forced to retreat all the way back to the first vertex...
            {
                cl = placesbeen.top();
                // cout << "From the top, our current location is " << cl << endl;
                for (int index1 = 0;
                     index1 < thegraph[cl].size(); index1++) // for each place we can go from our current location...
                {
                    prosp = thegraph[cl].at(index1);
                    //cout << "The prospective vertex we consider is " << prosp << endl;
                    if (beenthere[prosp] == 0) // if we haven't been there before...
                    {
                        // cl = prosp; // move our current location to there
                        //cout << "We are moving to " << prosp << endl;

                        placesbeen.push(prosp);      // and add it to the path
                        //cout << "Current location after pushing: " << placesbeen.top() << endl;
                        if (pc[placesbeen.top()] >
                            0) { return 1; }        // If our new site has a pebble, quit looking for pebbles
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
                {                            // Then we didn't move anywhere, so we need to pop off the last vertex and retreat
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
    bool SiteRP::loadsite(int i) {
        if (placesbeen.size() == 0) {
            while (pc[i] < 2 && findpebble(i)) // while the site is not loaded and you are finding pebbles...
                // c++ documentation says && short circuits, so you shouldn't
                // even look for a pebble if the site is loaded.
            {
                reversepath();                    // reverse the path, which also shifts the pebbles
            }
        }
        else {
            std::cout << "I tried to find a pebble, but when I started the path stack wasn't empty, it had size " <<
            placesbeen.size() << std::endl;
        }

        if (pc[i] == 2) { return 1; }        // return 1 if the site loaded successfully.
        else { return 0; }
    }

// loadsite with two arguments does the same thing, but won't try to take pebbles from skip to move onto i
// Which is to keep the two sites from swapping the three pebbles between themselves endlessly
    bool SiteRP::loadsite(int i, int skip) {
        if (placesbeen.size() == 0) {
            while (pc[i] < 2 && findpebble(i, skip)) // while the site is not loaded and you are finding pebbles...
                // C++ documentation says && short-circuits, so this shouldn't even
                // look for a pebble if the site is loaded.
            {
                reversepath();                    // reverse the path, which also shifts the pebbles
            }
        }
        else {
            std::cout << "I tried to load site " << i << " but the placesbeen stack wasn't empty.\n";
        }
        if (pc[i] == 2) { return 1; }        // return 1 if the site loaded successfully.
        else { return 0; }
    }

// loadsites tries to move pebbles until there are two on both sites i and j
    bool SiteRP::loadsites(int i, int j) {

        while (pc[j] < 2 && findpebble(j)) {
            reversepath();
        }

        while (pc[i] < 2 && findpebble(i, j))    // while the first site is not loaded and you are finding pebbles skipping
        {                                        // second site, load the first site, then load the second.
            reversepath();
            while (pc[j] < 2 && findpebble(j)) {
                reversepath();
            }
        }
        if (pc[i] == 2 && pc[j] == 2) { return 1; }
        else { return 0; }
    }

// addbond tries to load the sites. If it succeeds, it adds an edge from i to j and takes a pebble from i. Otherwise, it adds a redundant edge
    void SiteRP::addbond(int i, int j) {
        if (numbonds < 2 * size - 3 &&
            loadsites(i, j))            // If there are at least four pebbles left, we try to load the sites
        {
            addedge(i, j);                // if we succeed, we add the edge from i to j and remove a pebble from i
            numbonds++;
            edges.push_back(Bond(i, j));
            pc[i]--;
        }
        else {
            addredundant(i, j);            // otherwise, we leave the pebbles where we shuffled them and
            rbonds++;
        }                                // place only a redundant bond
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// REPORTING TO SCREEN

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void SiteRP::setfilestream(float cval, int tval) {


        char mychar[] = "./data/cxxxtxxxx.txt";
        int start = sizeof(mychar) - 13;
        cval = cval * 100;
        int intc = cval;
        mychar[start + 4] = '0' + (tval % 10000) / 1000;
        mychar[start + 5] = '0' + (tval % 1000) / 100;
        mychar[start + 6] = '0' + (tval % 100) / 10;
        mychar[start + 7] = '0' + (tval % 10);

        mychar[start] = '0' + (intc % 1000) / 100;
        mychar[start + 1] = '0' + (intc % 100) / 10;
        mychar[start + 2] = '0' + (intc % 10);


        char mychar_rcluster[] = "./data/rcluster_cxxxtxxxx.txt";
        int start_rcluster = sizeof(mychar_rcluster) - 13;
        mychar_rcluster[start_rcluster + 4] = '0' + (tval % 10000) / 1000;
        mychar_rcluster[start_rcluster + 5] = '0' + (tval % 1000) / 100;
        mychar_rcluster[start_rcluster + 6] = '0' + (tval % 100) / 10;
        mychar_rcluster[start_rcluster + 7] = '0' + (tval % 10);

        mychar_rcluster[start_rcluster] = '0' + (intc % 1000) / 100;
        mychar_rcluster[start_rcluster + 1] = '0' + (intc % 100) / 10;
        mychar_rcluster[start_rcluster + 2] = '0' + (intc % 10);

        myfile.close();
        myfile.open(mychar);

        rclusterfile.close();
        rclusterfile.open(mychar_rcluster);

        //cout << mychar;
    }

    void SiteRP::log(int span) {
        myfile << ll << "\t" << correlation << "\t" << numparts << "\t" << numbonds << "\t" << rbonds << "\t" <<
        giantsize_bond << "\t" << giantsize_site << "\t" << span << "\n";
    }

// listedges lists the edges from site i
    void SiteRP::listedges(int i) {
    std::cout << "\nvertex number " << i << " points towards the following vertices:\n";

        for (int index = 0; index < thegraph[i].size(); index++) {
            std::cout << thegraph[i].at(index) << " ";
        }

    std::cout << std::endl;
    }

// listalledges lists all the edges from the sites
    void SiteRP::listalledges() {
        for (int index = 0; index < size; index++) {
            listedges(index);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// FIXED-PROBABILITY TRIALS AND ATTENDANT FUNCTIONS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// BOND DILUTED NETWORK TO CONSIDER RIGIDITY PERCOLATION

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// functions in original codes by Zeb


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// SITE DILUTED NETWORK TO CONSIDER RIGIDITY PERCOLATION

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// initgiantrigidcluster() initializes the giantrigidcluster graph
    void SiteRP::initgiantrigidcluster() {
        giantsize_site = 0;
        giantsize_bond = 0;

        // rewrite the RigidIndex for edges
        for (std::vector<Bond>::iterator it = edges.begin(); it != edges.end(); ++it) {
            it->initBondRigidIndex();
        }

        for (int bondindex = 0; bondindex < size; bondindex++)  // Clear the graphs of giantrigidcluster
        {
            giantrigidcluster[bondindex].clear();
        }
    }

// initemptytrigraph() updates numbonds, rbonds, thegraph, rgraph, placesbeen to a triangular graph with no particles or bonds
    int SiteRP::initemptytrigraph() {
        numbonds = 0; // the number of bonds
        edges.clear();
        numparts = 0; // the number of particles
        rbonds = 0;   // the number of redundant bonds
        SpanLastStatus = 0; // the initial spanning status is NO

        for (int pcindex = 0; pcindex < size; pcindex++) // Just setting the pebble count to 2 everywhere.
        {                                               // and setting which sites are occupied
            pc[pcindex] = 2;
            occ[pcindex] = 0;
        }

        while (placesbeen.size() > 0)                    // Clear the places been stack
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

    void SiteRP::addtricluster2(int site, float c) // Has already added the rigidcluster function, as well as the spanning cluster
    {

        if (occ[site] == 0) {
            int newsite;
            int numneighbors = occ[choosedir(site, 1)] + occ[choosedir(site, 2)] + occ[choosedir(site, 3)] +
                               occ[choosedir(site, 4)] + occ[choosedir(site, 5)] + occ[choosedir(site, 6)];

            if (randprob() < pow(1. - c, 6 - numneighbors)) {
                occ[site] = 1;

                numparts++;
                for (int k = 1; k <= 6; k++) {
                    newsite = choosedir(site, k);
                    if (occ[newsite] == 1 && isempty(site, newsite)) {
                        addbond(site, newsite);
                    }
                }
                // choose some densities for the rigid cluster
                if (numparts % (ll / 2) == 0) {
                    rigidcluster();

                    int span = spanningrcluster();
                    log(span);

                    // choose the critical position to get rigid cluster decomposition
                    if (SpanLastStatus == 0 && span == 1){
                        StoreRigidInfoOfSite();
                    }

                    // update the SpanLastStatus
                    SpanLastStatus = span;
                }
                else
                    log();

            }
        }
    }

    void SiteRP::onetritrial2(long long int maxout, float c) {
        int numattempts = 0;
        initemptytrigraph();

        while (numattempts < maxout && numparts < ll * ll) {
            numattempts++;
            addtricluster2(randsite0(), c);
            //plot();
        }

        myfile.close();
        rclusterfile.close();
    }

    void SiteRP::multictrial(long long int maxout, float c1, float c2, float dc, int numtrials) {
        for (correlation = c1; correlation <= c2; correlation += dc) {
            for (int mtc = 14; mtc <= numtrials; mtc++) {
                int numattempts = 0;
                initemptytrigraph();

                setfilestream(correlation, mtc);

                while (numattempts < maxout && numparts < ll * ll) {
                    numattempts++;

                    addtricluster2(randsite0(), correlation);
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// rigid cluster

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool SiteRP::isredundant(int i, int j)  //see if the test bond between (i,j) is redundant(dependent)
    {
        if (numbonds < 2 * size - 3 &&
            loadsites(i, j))            // if there are at least four pebbles left, we try to load the sites
        {
            return 0;                // if we succeed, then the edge from i to j is independent
        }
        else {
            return 1;            // otherwise, the edge is redundant(dependent)
        }
    }

    bool SiteRP::isbondrigid(Bond &a, Bond &b) //see if the two bonds a and b are rigid to each other
    {
        if (isredundant(a.vertices.first, b.vertices.first) && isredundant(a.vertices.first, b.vertices.second) &&
            isredundant(a.vertices.second, b.vertices.first) && isredundant(a.vertices.second, b.vertices.second)) {
            return 1; // if all the four test bonds between a1,b1; a1,b2; a2,b1; a2,b2 are redundant
        }
        else {
            return 0;
        }
    }

    void SiteRP::rigidcluster() // mark the rigid clusters
    {
        initgiantrigidcluster(); //empty the vector array first

        int rcnum = 0; //index of the rigid cluster
        giantsize_bond = 1; //if there exists any bond, the smallest possible giant rigid cluster bond size is 1
        giantindex = 0; //the index for the giant rigid cluster (in this function)

        for (std::vector<Bond>::iterator refBond = edges.begin(); refBond != edges.end(); ++refBond) {
            if (refBond->RigidIndex == 0) {
                int rclustersize_bond = 1;//the bond size of this rigid cluster
                rcnum++;
                refBond->RigidIndex = rcnum;
                for (std::vector<Bond>::iterator testBond = edges.begin(); testBond != edges.end(); ++testBond) {
                    if (testBond->RigidIndex == 0) {
                        //the test and ref bonds are not in some rigid clusters
                        if (isbondrigid(*refBond, *testBond) == true) { //if test bond is rigid with respect to refbond
                            testBond->RigidIndex = rcnum;
                            rclustersize_bond++;
                        }
                    }
                }
                if (rclustersize_bond >= giantsize_bond) { //update the size of the giant rigid cluster
                    giantsize_bond = rclustersize_bond;
                    giantindex = rcnum;
                }
            }
        }

        // pick out the giant rigid cluster and store it in the vector "giantrigidcluster" (!!! we need it to become a undirected adjacent list)

        for (std::vector<Bond>::iterator it = edges.begin(); it != edges.end(); ++it) {
            if (it->RigidIndex == giantindex) {
                int site_I = it->vertices.first;
                int site_J = it->vertices.second; //the two sites of the rigid bond
                giantrigidcluster[site_I].push_back(site_J);
                giantrigidcluster[site_J].push_back(site_I);
            }
        }

        for (int i = 0; i <= size - 1; i++) {
            if (!giantrigidcluster[i].empty()) {
                giantsize_site++;
            }
        }
    }


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// THE RIGID CLUSTER INFO (FOR CLUSTER STATISTICS INFO)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // store the rigid cluster decomposition information in sites (initial rcluster_site before the running of this function)
    void SiteRP::StoreRigidInfoOfSite(){
        // clear the stored rigid cluster info, which is already not useful
        for (int i = 0; i <= size - 1; ++i) {
            rcluster_site[i].clear();
        }

        // store the rigid cluster info to rcluster_site
        for (std::vector<Bond>::iterator it = edges.begin(); it != edges.end() ; ++it) {
            int site_I = it->vertices.first;
            int site_J = it->vertices.second; //the two sites of the rigid bond

            // if rcluster_site has not stored the RigidIndex, store it
            if (find(rcluster_site[site_I].begin(),rcluster_site[site_I].end(),it->RigidIndex) == rcluster_site[site_I].end()){
                rcluster_site[site_I].push_back(it->RigidIndex);
            }
            if (find(rcluster_site[site_J].begin(),rcluster_site[site_J].end(),it->RigidIndex) == rcluster_site[site_J].end()){
                rcluster_site[site_J].push_back(it->RigidIndex);
            }
            //push the rigid indices to the rcluster_site in order to store the info
        }

        // print out the rigid cluster decomposition info
        for (int st = 0; st <= size - 1; ++st) {
            if (rcluster_site[st].empty()){
                rclusterfile << "0" << "\n";
            }
            else{
                for (std::vector<int>::iterator it = rcluster_site[st].begin(); it != rcluster_site[st].end(); ++it) {
                    rclusterfile << *it << "\t";
                }
                rclusterfile << "\n";
            }
        }
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// SPANNING RIGID CLUSTER

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// We need to pick out the giant rigid cluster from the network and then determine if it is the spanning cluster.

    bool SiteRP::spanningrcluster() {
        // Do the DFS in the giant rigid cluster
        int beenthere[size] = {};
        std::fill_n(beenthere, size, EMPTY);
        std::stack<int> DFS_rcluster;
        // Find a starting point
        int v_start;

        for (v_start = 0; v_start < size - 1; v_start++) {
            if (!giantrigidcluster[v_start].empty()) {
                break;
            }
        }

        //v_start is the starting point for the DFS
        beenthere[v_start] = 0;
        DFS_rcluster.push(v_start);
        int cl; //cl is the current location
        int prosp; //prosp is the prospective location we're going to move to
        while (DFS_rcluster.size() > 0) // make sure the DFS is in the giant rigid cluster
        {
            cl = DFS_rcluster.top();
            for (int index1 = 0; index1 < giantrigidcluster[cl].size(); index1++) {
                prosp = giantrigidcluster[cl].at(index1);
                if (beenthere[prosp] == EMPTY) {
                    DFS_rcluster.push(prosp);
                    //Do sth here
                    //beenthere[prosp]=1;


                    //Get the x-displacement for prosp from the x-displacement for cl (x_prosp=x_cl +1 OR +0 OR -1)
                    if (prosp == choosedir(cl, 2) || prosp == choosedir(cl, 5)) {
                        beenthere[prosp] = beenthere[cl];
                    } else if (prosp == choosedir(cl, 1) || prosp == choosedir(cl, 6)) {
                        beenthere[prosp] = beenthere[cl] + 1;
                    } else if (prosp == choosedir(cl, 3) || prosp == choosedir(cl, 4)) {
                        beenthere[prosp] = beenthere[cl] - 1;
                    }
                    else {
                        std::cout << "Whoooops, sth goes wrong." << std::endl;
                    }

                    break;
                }
                else //Existing marked displacement can be compared to current place to see if the spanning cluster exists.
                {
                    if (abs(beenthere[cl] - beenthere[prosp]) == ll - 1) {
                        return true;
                    }
                }
            }
            if (cl == DFS_rcluster.top() || DFS_rcluster.size() == 0) {
                DFS_rcluster.pop();
            }
        }
        return false;// After the DFS still no finding about the spanning cluster, then it is not spanning.
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// PLOT THE NETWORK CONFIGURATION WITHOUT ADDING RIGID CLUSTER INFORMATION

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void SiteRP::onetritrial2_plot(float p, float c) {
        initemptytrigraph();
        //while(numparts <= p*ll*ll)
        while (numparts <= p * ll * ll) {
            addtricluster2_withoutRIGID(randsite0(), c);
        }
        myfile.close();
        rclusterfile.close();
    }

    void SiteRP::addtricluster2_withoutRIGID(int site, float c) // Has not added the rigidcluster function, as well as the spanning cluster
    {

        if (occ[site] == 0) {
            int newsite;
            int numneighbors = occ[choosedir(site, 1)] + occ[choosedir(site, 2)] + occ[choosedir(site, 3)] +
                               occ[choosedir(site, 4)] + occ[choosedir(site, 5)] + occ[choosedir(site, 6)];

            if (randprob() < pow(1. - c, 6 - numneighbors)) {
                occ[site] = 1;

                numparts++;
                for (int k = 1; k <= 6; k++) {
                    newsite = choosedir(site, k);
                    if (occ[newsite] == 1 && isempty(site, newsite)) {
                        addbond(site, newsite);
                    }
                }
                //log();
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// CRITICAL EXPONENTS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// correlation length exponent \nu
// \nu is given from different critical point position to look at fluctuation

// cluster-statistics exponent \tau
// at the transition point, n_s(\phi_P) ~ s^{-\tau}, so to get \tau, need to statistics for cluster sizes at transition point

// fractal dimension d_f: for mass scaling analysis, average mass for rigid spanning cluster, M ~ L^{d_f-2}
// d_{BB}: for backbone dimension, average mass of stressed backbone, M' ~ L^{d_{BB}-2}
// for site percolation, only one exponent for the fractal dimension, d_{BB} (actually d_{BB} and d_f are describing the same properties)

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// UNI-TEST & APPLICATIONS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void SiteRP::OneTrialTest(float cfor, int trial) //Generate one-time trial for triangular lattice (site RP)
    {
        //cfor = 0.0; // correlation constant
        //trial = 1; // trial counting

        setfilestream(cfor,trial);
        onetritrial2(ll*ll*100000000L,cfor);
    }

    void SiteRP::MultiTrialTest() //Generate multiple-times trial for triangular lattice (site RP)
    {
        multictrial(ll*ll*100000000L,0.0,0.1,.2,20); //That's it!
    }

    void SiteRP::PlotNetworkTest() //Generate network plot file
    {
        float cfor;
        float p;
        std::cout<<"Type in the number for constant c: \n";
        std::cin>>cfor;
        std::cout<<"Type in network density p: \n";
        std::cin>>p;


        onetritrial2_plot(p,cfor);
        rigidcluster();

        int span = spanningrcluster();

        std::cout<<"The size of the giant rigid cluster is "<< giantsize_site <<" with the spanning status in "<<span<<"\n";

        std::ofstream outfile;
        outfile.open("rclusterout.txt");
        for(int i=0; i < size; i++)
        {
            if(occ[i] != 0)
                outfile << i%ll << "\t" << i/ll << "\t" << !giantrigidcluster[i].empty() << "\n";
        }
        outfile.close();
    }


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// END OF THE CLASS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
