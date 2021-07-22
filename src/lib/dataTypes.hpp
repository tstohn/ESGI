#pragma once

#include <iostream>
#include <unordered_set>
#include <string>
#include <string.h>
#include <cstdio>

/**
 *    small functionality of a char* set to not use any duplicate strings
 *    
 * */

class CharHash
{
    public:
        size_t operator()(const char *s) const
        {

            return std::hash<std::string_view>()(std::string_view(s, std::strlen(s)));
        }
};

struct CharPtrComparator {
   bool operator()(const char* left, const char* right) const {
      //return ((left != nullptr) && (right != nullptr) && (strcmp(left, right) < 0));
      return ( (strcmp(left, right)) == 0 );

   }
};

class UniqueCharSet
{

   public:

      ~UniqueCharSet()
      {
         clearUniqueSet();
      }

      void printSet() 
      {
         for(auto it = charPtrSet.begin(); it != charPtrSet.end(); ++it)
            std::cout << *it << "\n";
      }

      const char* getUniqueChar(const char* k)
      {
         std::unordered_set<const char*>::iterator idx = charPtrSet.find(k);
         if(idx != charPtrSet.end())
         {
            return *idx;
         }
         else
         {
            const char* x = insertElement(k);
            return(x);
         }
      }

      const char* insertElement(const char* k) 
      {
         if(!k) {
            std::cout << "Unable to insert char*: \'" << k << "\'\n";
            exit(EXIT_FAILURE);
         }
         char* key = new char[strlen(k) + 1];

         std::cout << "ADRESS: " << static_cast<void*>(key) << " ";
         strcpy(key, k);

         std::cout << "\nITER THROUGH <" << *key << ">:_";
         for(char* x = key; *x !='\0'; ++x)
         {
            std::cout << *x;
         }
                  std::cout << "_\n";

         std::cout << "ADRESS: " << static_cast<void*>(key) << "\n";

         std::cout << "make new " << charPtrSet.size() << " } <"<<  k << ">|<" << key << ">\n";

         auto x = charPtrSet.insert(key);
         int o = 0;
         std::cout << "size: " << charPtrSet.size() << " " << *(x.first)  << " " << x.second<< true <<"\n";
         if(x.second == false)
         {
            std::cout << "EXIT" << "\n";
            exit(0);
         }
         //for(auto el : charPtrSet)
         //{
           // std::cout << o<< " "<< el << "\n";
          //  ++o;
         //}
         std::cout << "make new-@@@\n";

         return(key);
      }

      void clearUniqueSet() 
      {
         std::cout << "CLEARING SET\n";
         for(auto it = charPtrSet.begin(); it != charPtrSet.end(); ++it) 
         {
            const char* key = *it;

            if(key) {
               delete [] key;
               key = nullptr;
            }
         }
         charPtrSet.clear();
      }

   private:
      std::unordered_set<const char*, CharHash, CharPtrComparator> charPtrSet;

};