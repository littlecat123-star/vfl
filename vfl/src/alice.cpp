//This realizes the paper's online phase
#include "ttp.h"
#include "hash.h"
#include "inv.h"

#include "PackedArray.h"  // 
#include <omp.h>

using namespace emp;
using namespace std;
int port, party;

int main(int argc, char** argv) {
  omp_set_num_threads(THREADS);
  
    if (argc!=3) {
        cout <<"Must supply party and port number"<<endl;
    }

    parse_party_and_port(argv, &party, &port);

    cout <<"Buckets: "<<BUCKETS<<", BETA: "<<BETA<<", #Hash elements: "<<HASHELEMENTS<<", modulus: "<<MODULUS<<endl;
    
    
    //Only known to Alice
    BT *sA = (BT *) calloc(1, BUCKETS*sizeof(BT));
    BT *rA = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT));
    PackedArray* toBobPacked = NULL;  
    PackedArray* toAndyPacked = NULL; // NEW,send the finnal payload
    BT *aliceCuckooTable = NULL;
    BT *toBob = NULL;
    BT *fromBob = NULL;
    PackedArray* fromBobPacked = NULL;
    //new
    BT *fromBobPayload = NULL;
    PackedArray* fromBobPayloadPacked = NULL;
    BT *aliceTable2 = NULL; //save the finall payload table
    // before two item are used to save the payload from bob

    //let BOb be the corrordiantor (NEW)
    BT *bobTable2 = NULL;//save the 2-th ID hash
    // end new

    //Only known to Bob 
    BT *INV = (BT*) calloc(1, MODULUS*sizeof(BT));
    BT *bobTable = NULL; // save the 1-th ID hash
    BT *rB = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT));
    BT *sB = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT));
    PackedArray* fromAlicePacked = NULL;
    
    BT *fromAlice = NULL;
    PackedArray* toAlicePacked = NULL;
    BT *toAlicePayload = NULL;
    PackedArray* toAlicePayloadPacked = NULL; //save corresponding payload

    auto quad = clock_start();
    if (party == BOB) {
      computeInverses(INV);
      cout <<"Computing inverses in "<<time_from(quad)/1000<<" ms"<<endl;
    }
    //For simplicity, we quickly generate the same triples on the fly on both parties. This is typically done offline before.
    cout <<"Computing triples "<<flush;
    quad = clock_start();
    computeTriples(rA, sA, rB, sB);
    cout <<"in "<<time_from(quad)/1000<<" ms"<<endl;

    if (party == ALICE) {
      free(rB);
      free(sB);
      free(INV);
    } else {
      free(rA);
      free(sA);
      }
       
    //Fix seed for debugging of fixed input
    block seed = makeBlock(7, 8);
    PRG prg;
    prg.reseed(&seed);

    //We use Cuckoo hashing with 3 hashfunctions. We can prepare these offline.
    BT *h1 = (BT *) calloc(1, HASHELEMENTS*sizeof(BT));
    BT *h2 = (BT *) calloc(1, HASHELEMENTS*sizeof(BT));
    BT *h3 = (BT *) calloc(1, HASHELEMENTS*sizeof(BT));
    
    cout <<"Preparing hash functions..."<<endl;
    auto hashFTime = clock_start();
    computeHashFunctions(prg, h1, h2, h3);
    cout <<"Time to prepare hash functions: "<<time_from(hashFTime)/1000<<" ms"<<endl;
    
    double hashing = 0.0;
    double CPUTime = 0.0;
    
    if(party == ALICE) {
      //Now we start preparing the Cuckoo table
        toBobPacked = PackedArray_create(BITS, BUCKETS);
        aliceCuckooTable = (BT *) calloc(1, BUCKETS*sizeof(BT));
        toBob = (BT *) calloc(1, BUCKETS*sizeof(BT));
        BT *aliceInput  = (BT*) calloc(1, NELEMENTS * sizeof(BT));
        fromBob = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT));
        fromBobPacked = PackedArray_create(BITS, BETA*BUCKETS);
        fromBobPayload = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT)); // NEW
        fromBobPayloadPacked = PackedArray_create(BITS, BETA*BUCKETS); // NEW
        aliceTable2 = (BT *) calloc(1, 3*BUCKETS*sizeof(BT)); // new finally table size,3 hash functions
        //Fill Alice with fixed input
	for (size_t i = 0; i<NELEMENTS; i++) {
	  aliceInput[i] = i;
	}
        cout<<"NELEMENTS="<<NELEMENTS<<endl;
	
        cuckooEntry **bufAliceCuckooTable = (cuckooEntry **) calloc(1, BUCKETS*sizeof(cuckooEntry*));

        auto start = clock_start();

	//Populate Cuckoo table...
        computeCuckoo(aliceInput, h1, h2, h3, bufAliceCuckooTable);

        //...and convert it into a proper format to send it
        for (size_t i = 0; i < BUCKETS; i++) {
            if (bufAliceCuckooTable[i]!=NULL) {
                aliceCuckooTable[i] = (bufAliceCuckooTable[i])->value;
            } else {
	      //Mark as invalid
                aliceCuckooTable[i] = 3*HASHELEMENTS;
            }
        }

	hashing = time_from(start);
	CPUTime += hashing;
        cout <<"Alice Cuckoo hashing time: "<<hashing/1000<<" ms"<<endl;

	//No use for hash functions anymore
	free(h1);
	free(h2);
	free(h3);

	free(bufAliceCuckooTable);
	
        //Compute $c$s for Bob
	double elapsed = 0.0;
        start=clock_start();
        for (size_t i = 0; i < BUCKETS; i++) {
            toBob[i] = (sA[i] - aliceCuckooTable[i] + MODULUS) % MODULUS;
        }
        PackedArray_pack(toBobPacked, 0, toBob, BUCKETS);
	elapsed = time_from(start);
	CPUTime += elapsed;
        cout <<"Alice sA - H(x) time: "<<elapsed/(1000)<<" ms"<<endl;

	free(aliceCuckooTable);
	free(aliceInput);
	
        //END ALICE
    } else {
        //BOB

      //We prepare a regular hash table
        bobTable = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT));
        // bobTable2 = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT)); // init the second hash table
        fromAlicePacked = PackedArray_create(BITS, BUCKETS);
        PackedArray* fromAlicePayloadPacked = NULL; // new as Andy
        fromAlicePayloadPacked = PackedArray_create(BITS, 3*BUCKETS);//new
        fromAlice = (BT *) calloc(1, BUCKETS*sizeof(BT));
        toAlicePacked = PackedArray_create(BITS, BETA*BUCKETS);
        toAlicePayload = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT)); // new
	toAlicePayloadPacked = PackedArray_create(BITS, BETA*BUCKETS); // save the payload
        //Compute fixed input
        BT *bobInput  = (BT*) malloc(2 * NELEMENTS * sizeof(BT));
        for (size_t i = 0; i<2* NELEMENTS; i++) {
	  bobInput[i] = i;
	}
	cout<<"NELEMENTS="<<NELEMENTS<<endl;
        bucket *_bobTable = (bucket*) calloc(1, BUCKETS*sizeof(bucket));
        bucket *_toAlicePayload = (bucket*) calloc(1, BUCKETS*sizeof(bucket)); // new
        // bucket *_bobTable2 = (bucket*) calloc(1, BUCKETS*sizeof(bucket)); // new
        auto t = clock_start();

	//Compute regular hash table
	computeHash(bobInput, h1, h2, h3, _bobTable, bobTable,_toAlicePayload, toAlicePayload); // bobTable is the simple hash table
	// computeHash(bobInput, h1, h2, h3, _bobTable2, bobTable2,); // new
	hashing = time_from(t);
	CPUTime += hashing;
        cout <<"Bob hashing and computing payload: "<<hashing/1000<<" ms"<<endl; // double time
    
	//No use for hash functions anymore
	free(h1);
	free(h2);
	free(h3);

	free(_bobTable);
        free(_toAlicePayload);
        // free(_bobTable2);
	free(bobInput);
        //END BOB
    }

    cout <<"Connecting..."<<endl;
    double total_time = 0.0;
    auto ttime_including_idle = clock_start();  // waiting time

    //Server or local environment: NetIO* io = new NetIO(party==ALICE?nullptr:"127.0.0.1", port);
    //West Coast: NetIO* io = new NetIO(party==ALICE?nullptr:"104.42.191.30", port);
    //EU: NetIO* io = new NetIO(party==ALICE?nullptr:"13.95.162.209", port);

    
    NetIO* io = new NetIO(party==ALICE?nullptr:"127.0.0.1", port);
    
    if (party == ALICE) {

      //Send the $c$s
      
        //4 = number of bytes per uint32_t slot
        size_t myLength = 4*PackedArray_bufferSize(toBobPacked);
        size_t myLength2 = 4*PackedArray_bufferSize(toBobPacked); // new
        cout <<"Alice sending "<<myLength<<" Bytes"<<endl;
        io->sync();
	
	ttime_including_idle = clock_start();
        auto start=clock_start();
        mySend(io, (char*)toBobPacked, myLength);
        auto elapsed = time_from(start);
        total_time +=elapsed;
        cout <<"Alice sent c in "<<elapsed/1000<<" ms"<<endl;

	free(toBobPacked);
        
	//Alice gets back $d$s
        myLength = 4*PackedArray_bufferSize(fromBobPacked);
        myLength2 = 4*PackedArray_bufferSize(fromBobPayloadPacked);
        
        cout <<"Alice receiving "<<myLength+myLength2<<" Bytes"<<endl; // double length
        io->sync();     
        start=clock_start();
        myRecv(io, (char*) fromBobPacked, myLength);       
	auto pack = clock_start();
        PackedArray_unpack(fromBobPacked, 0, fromBob, BETA * BUCKETS);
        CPUTime += time_from(pack);	
        elapsed = time_from(start);
        total_time += elapsed;
        cout<<"Alice receive d"<<endl;
        cout <<"Alice received d in "<<elapsed/1000<<" ms"<<endl;
        io->sync();     
        start=clock_start();
        myRecv(io, (char*) fromBobPayloadPacked, myLength2); //new
        auto packpayload = clock_start();
        PackedArray_unpack(fromBobPayloadPacked, 0, fromBobPayload, BETA * BUCKETS);
	CPUTime += time_from(packpayload);	
        elapsed = time_from(start);
        total_time +=elapsed;
        cout<<"Alice receive payload"<<endl;
        cout <<"Alice received payload in "<<elapsed/1000<<" ms"<<endl;
	free(fromBobPacked);
        free(fromBobPayloadPacked); // new
        //Compute intersection
        start=clock_start();
        size_t match = 0;
        cout<<"Alice compute psi"<<endl;
	#pragma omp parallel for reduction(+:match)
        for (size_t i = 0; i < BUCKETS; i++) {
	  for (size_t j = 0; j < BETA; j++) {
                if (rA[i*BETA + j]==fromBob[i * BETA + j]) {
                    match++;
                    // new
                    for (size_t k = 0; k < 3; k++) 
                    {
		       aliceTable2[i*3 + k]=fromBobPayload[i*BETA + j]+k ; 
                    }            
                    }
            }
        }
        toAndyPacked = PackedArray_create(BITS, BUCKETS*3);
        PackedArray_pack(toAndyPacked, 0, aliceTable2, BUCKETS);
        elapsed = time_from(start);
        total_time +=elapsed;
	CPUTime+=elapsed;
        cout <<"Alice computes the A-C psi and payload in "<<elapsed/1000<<" ms, size "<<match<<endl;
	// send the finnal payload to compare  
        myLength = 4*PackedArray_bufferSize(toAndyPacked); 
        cout <<"Alice sending the finnal payload "<<myLength<<" Bytes"<<endl;      
        io->sync();
        start=clock_start();
	mySend(io, (char*) toAndyPacked, myLength);	
        elapsed = time_from(start);
        total_time +=elapsed;
        cout <<"Alice sent finall payload in "<<elapsed/1000<<" ms"<<endl;
	free(toAndyPacked);
	//END ALICE
    } else {
        //BOB

      //Bob gets the $c$s
      
        //4 = number of bytes per uint32_t slot
        size_t myLength = 4*PackedArray_bufferSize(fromAlicePacked);
        size_t myLength2 = myLength; // new 
        //BOB
        cout <<"Bob receiving "<<myLength<<" Bytes"<<endl;
        io->sync();
	
	ttime_including_idle = clock_start();
        auto start=clock_start();

        //io->recv_data(fromAlicePacked, myLength);
	myRecv(io, (char *) fromAlicePacked, myLength); // cannot receive?
        cout<<"receiving c"<<endl;
	auto unpack = clock_start();
        PackedArray_unpack(fromAlicePacked, 0, fromAlice, BUCKETS);
        cout<<"unpack the c"<<endl;
	auto elapsed = time_from(start);
	CPUTime +=time_from(unpack);
        total_time +=elapsed;
        cout <<"Bob received in "<<elapsed/1000<<" ms"<<endl;

	free(fromAlicePacked);
        
        //Compute result: d =  (fromAlice + myHash + sB) / rB 
        start = clock_start();
	bobsOperations(bobTable, fromAlice, rB, sB, INV); // compute d,d is saved in bobTable
	PackedArray_pack(toAlicePacked, 0, bobTable, BETA*BUCKETS); // pack
        cout<<"BOb pack d"<<endl;
	PackedArray_pack(toAlicePayloadPacked, 0, toAlicePayload, BETA*BUCKETS); // new
        cout<<"BOb pack payload"<<endl;
        elapsed = time_from(start);
        total_time +=elapsed;
	CPUTime+=elapsed;
        cout <<"Bob computes answer in "<<elapsed/1000<<" ms"<<endl;

	free(fromAlice);
	free(bobTable);
        free(toAlicePayload); // NEW


	//Send all $d$s back
        myLength = 4*PackedArray_bufferSize(toAlicePacked);
        myLength2 = 4*PackedArray_bufferSize(toAlicePayloadPacked); // new
        cout <<"Bob sending "<<myLength+myLength2<<" Bytes"<<endl;
        io->sync();
        start=clock_start();
	mySend(io, (char*) toAlicePacked, myLength);	
        elapsed = time_from(start);
        total_time +=elapsed;
        cout <<"Bob sent d in "<<elapsed/1000<<" ms"<<endl;
	free(toAlicePacked);
        io->sync();
        start=clock_start();
	mySend(io, (char*) toAlicePayloadPacked, myLength2);	
        elapsed = time_from(start);
        total_time +=elapsed;
        cout <<"Bob sent payload in "<<elapsed/1000<<" ms"<<endl;       
        free(toAlicePayloadPacked);


        // serve as Andy
        start=clock_start();
        PackedArray* fromAlicePayloadPacked = NULL; // new as Andy
        fromAlicePayloadPacked = PackedArray_create(BITS, 3*BUCKETS);//new    
        myLength = 4*PackedArray_bufferSize(fromAlicePayloadPacked);
        cout <<"Andy receiving "<<myLength<<" Bytes"<<endl;
        io->sync();  
        //io->recv_data(fromAlicePacked, myLength);
	myRecv(io, (char *) fromAlicePayloadPacked, myLength); // cannot receive?
        cout<<"Andy receiving payload from alice "<<endl;
        free(fromAlicePayloadPacked);
        elapsed = time_from(start);
        cout<<"Andy receive final payloadk in "<<elapsed/1000<<" ms"<<endl;


    }

    cout <<"Total time for party "<<myName<<": "<<total_time/1000<<" ms, sent: "<<(io->counter)<<" Byte, i.e., "<<((double)io->counter)/(1024*1024)<<" MByte"<<endl;
    // running time after hashing without waiting
    cout <<"Total time including waiting for "<<myName<<": "<<time_from(ttime_including_idle)/1000<<" ms"<<endl;
    // running time after hashing with waiting
    cout <<"CPU time for "<<myName<<": "<<(CPUTime)/1000<<" ms"<<endl;
    // computing time 
    // cout<<"hashing time "<<hashing/1000<<"ms"<<endl;
    cout <<"Total time including waiting and hashing for "<<myName<<": "<<(time_from(ttime_including_idle)+hashing)/1000<<" ms"<<endl;
    // hash time+running time
    
}
