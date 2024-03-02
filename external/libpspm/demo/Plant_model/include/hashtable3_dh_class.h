#ifndef HASHTABLE_CPU_GPU
#define HASHTABLE_CPU_GPU
#include <iostream>
#include <cassert>
/*using namespace std;*/

//#include "../utils/simple_math.h"

#define __DEBUG_HT_  if(false) 


//template <typename Key>
//struct KeyHasher{
//	int operator()(const Key key) const {
//		return (key.z*377771 + key.y*677 + key.x);	// 2.15, 64
//	}
//};


// HASH TABLE IMPLEMENTATION 

template <typename Key, typename Value>
class HashNode{
	public:
	Key key;
	Value value;			
	int count;				// Max number of probes ever required to store a key that hashed to this slot (the hashtable is just an array of HashNodes, so node=slot)
	bool isEmpty;			// Is this slot empty?

	bool isDeleted;			// for visualization purpose only.
	
	HashNode(){
		key = Key();		// Value-initialization for Key and Value. This will be 0 for built-in types and constructor-initialized for user defined types
		value = Value();
		count = 0;			
		isEmpty = true;		
		isDeleted = false;	
	}

	void print(){
		if (isDeleted) std::cout << "---<" << key << "," << value << ">--- " << " [" << count << " probes]" << std::endl;
		else if (!isEmpty) std::cout << "<" << key << ", " << value << "> " << " [" << count << " probes]" << std::endl;
		else std::cout << "-----------" << " [" << count << " probes]" << std::endl;
	}
};


template <typename Key>
int keyHasher(Key key, int length){
//	return ((unsigned long long)(log2(key+1))*377771) % length;
	return ((unsigned long long)(key)*377771) % length;
//	return ((unsigned long long)key*677) % length;
//	return (key.z*377771 + key.y*677 + key.x) % length;	// 2.15, 64
//	return (key.z*377771 + key.y*133337 + key.x) % length; // 1.7, 25
//	return (key.z*497771 + key.y*133337 + key.x) % length; // 30, 173
//	return ((long long int)key.z*497771 + key.y*787 + key.x) % length; // 1.5, 28
}


template <typename Key>
int incrementHasher(Key key){
//	return 1;
	return 2*((key*677)%53)+17;	// should return an odd number, as table size is 2^m
}


//template <typename Key>
//int increment(Key key, int count){
////	return (id+11)%length;
//	return count*incrementHasher(key);
////	return (key)%677;
//}

struct HashFindResult{
	size_t id;			// id where key was found (or -1 if not found)
	int    attempts;	// how many probes were required to find the key or conclude not-found?
};


template <typename Key, typename Value>
class HashTable{
	public: 
	int length;
	float load_factor;
	HashNode<Key,Value> * ht;
	int nelements = 0;
	
	HashTable(){
		ht = nullptr;
		__DEBUG_HT_ std::cout << "Constructor called @" << ht << std::endl;
	}
	
	HashTable(int len, float lf){
		assert (len > 0);
		length = len;
		load_factor = lf;
		ht = new HashNode<Key,Value> [len];
		__DEBUG_HT_ std::cout << "Constructor called @" << ht << std::endl;
	}	


	HashTable(HashTable& _h){
		length = _h.length;
		load_factor = _h.load_factor;
		ht = new HashNode<Key,Value> [length];
		for (int i=0; i<length; ++i) ht[i] = _h.ht[i];
	}


	HashTable& operator =(const HashTable& _h){
		length = _h.length;
		load_factor = _h.load_factor;
		ht = new HashNode<Key,Value> [length];
		for (int i=0; i<length; ++i) ht[i] = _h.ht[i];
		__DEBUG_HT_ std::cout << "Assignment called @" << ht << std::endl;
		return *this;
	}

	
	~HashTable(){
		__DEBUG_HT_ std::cout << "Destructor called @" << ht << std::endl;
		delete [] ht;
		ht = nullptr;
	}
	
	
	HashFindResult hash_find(Key key) const{
		HashFindResult res;
		
		size_t hash = keyHasher(key, length);
		size_t id;
		int count = 0; 
		__DEBUG_HT_ std::cout << "@" << ht << " > Find: (" << key << "): |";// << "*, ";
		while(1) { 	// probe until slot is filled AND the key matches. (i.e., dont stop at empty slot)
			id = (hash + count*incrementHasher(key))%length; 
			++count;


			__DEBUG_HT_{
				std::cout << id; 
				std::cout.flush();
				// if (id == res.firstEmptySlot) std::cout << "__";
				std::cout << ", ";
			}
			
			if ((ht[id].key == key && !ht[id].isEmpty)) break;	// loop breaks as soon as key is found AND the slot was not marked empty
			
		 	// if probe reaches this point on the count'th try, the key cannot be in the table, because it it were, loop would've exited on the previous statement. (Only 'count' no. of entries were ever stored for a key hasshing to here)
			if (count >= ht[hash].count) { id = -1; break; }	// this can be split into 2 statements, 1 to set id, 2nd break condition. 
		
			assert(count <= length);	// if count exeeds length, something is wrong! -- This CAN occur is probes accumulate, but rebuilding the HT will solve the problem
		} 

		res.attempts = count;
		res.id = id;
		
		__DEBUG_HT_ {
			if (id != size_t(-1)) std::cout << "|   <" << key << ", " << ht[id].value << "> [" << res.attempts << " tries]" << /*((key != ht[id].value*10)? " * FAIL *":"") <<*/ std::endl;
			else std::cout << "|   <" << key << ", " << "- " << "> [" << res.attempts << " tries]" << std::endl;
		}
		
		return res;
	}


	// Insertions will happen on CPU only (Dont want headache of thread locks)
	int hash_insert(Key key, Value value, int* duplicate=NULL){
		__DEBUG_HT_ std::cout << "Insert: <" << key << ", " << value << "> : trying "; std::cout.flush();
		size_t hash = keyHasher(key, length);	// get the hash of the specifed key

		size_t id;								// the location where this <key, value> will be stored equals hash, but will be incremented if not empty 
		size_t firstEmptySlot = -1;
		int    count_firstEmptySlot = length+1;	// set to index beyond the table's range
		int count = 0;							// number of collisions before empty slot was found (defaults to 0, i.e. no increments were required)
		do { 	
			id = (hash + count*incrementHasher(key))%length;			// increment slot (using double hashing) until empty or deleted slot found 
			++count;							// counter is not incremented if inserting into previously deleted slot, because keys that have been mapped past the deleted slot stay intact

			if (ht[id].isEmpty && firstEmptySlot == -1){	// if found an empty slot and firstEmptySlot is still invalid
				firstEmptySlot = id;
				count_firstEmptySlot = count;
			}
			
			__DEBUG_HT_{
				std::cout << id;
				if (id == firstEmptySlot) std::cout << "__";
				std::cout << ", ";
				std::cout.flush();
			}

		 	// if probe exceeds count, the key cannot be in the table, and we have reached a slot beyond count
			if (count > ht[hash].count) break;
			
			assert(count <= length);
			
		} while (!(ht[id].key == key && !ht[id].isEmpty));

		// If key was found in HT, we simply need to replace the associated value and exit
		if (ht[id].key == key && !ht[id].isEmpty){	
			__DEBUG_HT_ std::cout << "(count = " << count << ") --> Replacing @ [" << id << "]: <" << key << "," << ht[id].value << "> ==> <" << key << "," << value << ">" << std::endl;
			ht[id].value = value; 
			if (duplicate != NULL) *duplicate = 1;	// duplicate is used as a duplicate detection flag
			return count;
		}

		__DEBUG_HT_ std::cout << "... ";
		// Otherwise (if key was not found) and no empty slot encountered so far (firstEmptySlot != -1)
		//     continue searching for an empty slot 
		while (!ht[id].isEmpty && firstEmptySlot == -1){	
			assert(count <= length);
			id = (hash + count*incrementHasher(key))%length;			// increment slot (using double hashing) until empty or deleted slot found 
			++count;										// in this loop, counter must be incremented AFTER testing the ID, because previous loop exited after incrementing the counter. counter is not incremented if inserting into previously deleted slot, because keys that have been mapped past the deleted slot stay intact
			__DEBUG_HT_ std::cout << id << ", ";
		}
		__DEBUG_HT_ std::cout << "(count = " << count << ")";
		
		// If empty slot was detected, set ID and count to that of the empty slot 
		if (firstEmptySlot != -1){
			id = firstEmptySlot;
			count = count_firstEmptySlot;
		}
		
		ht[id].value = value;				// once empty slot found, store the <key,value> in the slot
		ht[id].key = key;
		ht[id].isEmpty = false;
		ht[id].isDeleted = false;
		ht[hash].count = std::max(ht[hash].count, count);	// when insertion happens in a deleted slot, we dont want to lose count of how many valid entries existed beyond the deleted slot. Hence, count must be maximum of previous count and result of current insertion
		__DEBUG_HT_ std::cout << " --> stored @ [" << id << "]" << std::endl;
		if (duplicate != NULL) *duplicate = 0;	// duplicate is used as a duplicate detection flag
	//	if (duplicate != NULL) *duplicate = id;

		++nelements;	// dont increment if (duplicate) key is replaced
		return count;

	}


	int hash_delete(Key key){
		size_t hash = keyHasher(key, length);
		auto res = hash_find(key);
		size_t id = res.id;
		int count = res.attempts;
		if (id == size_t(-1)) return 1;	// if key is not found, nothing to be done
		ht[id].isEmpty = true;
		ht[id].isDeleted = true;
		
		__DEBUG_HT_ std::cout << "Delete: <" << key << "," << ht[id].value << "> (count, attempts-1) = " << ht[hash].count << ", " << count-1 << ")" << std::endl;	// TODO: Can reduce count by one when the last number is deleted, but that may well have been the only number inserted, in which case count should have gone to zero. But theres no way to keep track of this! When creating the hashmap for interpolator, upon interval refinement, dont delete the interval. Only change the value, and add a new interval. 
		return 0;
	}


	void hash_clear(){
		for (int i=0; i<length; ++i) ht[i] = HashNode<Key,Value>();
	}


	void hash_print(){
		std::cout << "~~~ HASH TABLE ~~~~~~~~~~~" << std::endl;
		for (int i=0; i< std::min(length,10); ++i){
			std::cout << i << ": "; ht[i].print();
		}	
		std::cout << "Location: " << ht << std::endl;
		std::cout << "# Elements: " << nelements << std::endl;
		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	}

	
	
};









#endif






