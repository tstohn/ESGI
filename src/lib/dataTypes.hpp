#pragma once

#include <iostream>
#include <unordered_set>
#include <string>
#include <string.h>
#include <cstdio>

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
         std::unordered_set<const char*, CharHash, CharPtrComparator>::iterator idx = charPtrSet.find(k);
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

      void clearUniqueSet() 
      {
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

      const char* insertElement(const char* k) 
      {
         if(!k) {
            exit(EXIT_FAILURE);
         }
         char* key = new char[strlen(k) + 1];
         strcpy(key, k);
         charPtrSet.insert(key);

         return(key);
      }

      std::unordered_set<const char*, CharHash, CharPtrComparator> charPtrSet;
};