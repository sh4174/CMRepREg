/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: Registry.h,v $
  Language:  C++
  Date:      $Date: 2008/04/03 11:30:33 $
  Version:   $Revision: 1.1 $
  Copyright (c) 2003 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#ifndef __Registry_h_
#define __Registry_h_

#ifdef _MSC_VER
# pragma warning(disable:4786)  // '255' characters in the debug information
#endif //_MSC_VER

#include <cstdio>
#include <cstring>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <vector>

template <class T>
inline T GetValueWithDefault(const std::string &source, bool isNull, T defaultValue)
{
  // Initialize with the default value
  T returnValue = defaultValue;

  // Default value is returned if the entry is Null
  if(isNull)
    return returnValue;

  // Try to access the value using c++ formatting
  std::istringstream iss(source);
  iss >> returnValue;

  // Proceed as if the operation succeeded
  return returnValue;
}

template <>
inline const char *GetValueWithDefault<const char *>(const std::string &source, bool isNull, const char *defaultValue)
{
  if(isNull)
    return defaultValue;
  else
    return source.c_str();
}

/** A class used to associate a set of strings with an ENUM so that the
 * enum can be exported to the registry in a consistent fashion */
template <class TEnum>
class RegistryEnumMap
{
public:
  typedef std::string StringType;
  void AddPair(TEnum value, const char *description)
  {
    m_EnumToStringMap[value] = std::string(description);
    m_StringToEnumMap[description] = value;
  }
  
  const StringType operator() (TEnum eval) 
    { return m_EnumToStringMap[eval]; }

private:
  std::map<TEnum, StringType> m_EnumToStringMap;
  std::map<StringType, TEnum> m_StringToEnumMap;
  friend class RegistryValue;
};

/** A class that represents a value in the registry */
class RegistryValue
{
public:
  /** Default constructor: sets this object to Null state */
  RegistryValue();

  /** Initializing constructor: sets the object to a value */
  RegistryValue(const std::string &value);

  /** Is the value present in the Registry? */
  bool IsNull() const { return m_Null; }

  /** Get the internally stored string */
  const std::string &GetInternalString() const { return m_String; }

  /**
   * An operator that allows us to access a value using different
   * formats
   */
  bool operator[](bool defaultValue) {
    return GetValueWithDefault(m_String,m_Null,defaultValue);
  }

  int operator[](int defaultValue) {
    return GetValueWithDefault(m_String,m_Null,defaultValue);
  }

  unsigned int operator[](unsigned int defaultValue) {
    return GetValueWithDefault(m_String,m_Null,defaultValue);
  }

  double operator[](double defaultValue) {
    return GetValueWithDefault(m_String,m_Null,defaultValue);
  }

  const char *operator[](const char *defaultValue) {
    return GetValueWithDefault(m_String,m_Null,defaultValue);
  }

  std::string operator[](const std::string &defaultValue) {
    if(IsNull())
      return defaultValue;
    else
      return m_String;
  }

  /**
   * An operator that allows us to write a value using different formats
   */
  template <class T> void operator << (const T newValue)
  {
    // Create an output stream
    std::ostringstream oss;

    // Put the new value into the stream
    oss.precision(16);
    oss << newValue;

    // Set the string
    m_String = oss.str();
    m_Null = false;
  }

  /** Put an enum into this entry */
  template <class TEnum> void PutEnum(
    RegistryEnumMap<TEnum> &rem,TEnum value)
  {
    (*this) << rem.m_EnumToStringMap[value];;
  }

  /** Get an enum from this entry */
  template <class TEnum> TEnum GetEnum(
    RegistryEnumMap<TEnum> &rem,TEnum defaultValue)
  {
    if(rem.m_StringToEnumMap.find(m_String) == rem.m_StringToEnumMap.end())
      return defaultValue;
    else return rem.m_StringToEnumMap[m_String];
  }

private:
  std::string m_String;
  bool m_Null;
};




/**
 * \class Registry
 * \brief A tree of key-value pair maps
 */
class Registry
{
public:
  // String definition
  typedef std::string StringType;

  // List of strings
  typedef std::list<StringType> StringListType;

  /** Constructor initializes an empty registry */
  Registry();

  /** Constructor loads a registry from a file */
  Registry(const char *fname);

  /** Destructor */
  virtual ~Registry();

  /** Get a reference to a value in this registry, which can then be queried */
  RegistryValue &operator[](const StringType &key) { return Entry(key); }

  /** Get a reference to a folder inside this registry, creating it if necessary */
  RegistryValue &Entry(const StringType &key);

  /** Get a reference to a folder inside this registry, creating it if necessary */
  Registry &Folder(const StringType &key);

  /** Check whether there is a folder by a name or not */
  bool IsFolder(const StringType &key);

  /** A helper method to convert a printf-style expression to a key */
  static StringType Key(const char *key,...);

  /** Get a list of all keys that have values, append it to keyList */
  int GetEntryKeys(StringListType &keyList) const;

  /** Get a list of all subfolder keys, append it to keyList */
  int GetFolderKeys(StringListType &keyList) const;

  /** Find a value in a folder or return "" */
  StringType FindValue(const StringType& value);

  /** Empty the contents of the registry */
  void Clear();

  /** Get a list of all keys that have values contained in this registry
   * and all subfolders (recursive operation).  Prefix is used internally,
   * but can be specified to prepend a string to all keys */
  void CollectKeys(StringListType &keyList,const StringType &keyPrefix="");

  /** Remove all keys from this folder that start with a string */
  void RemoveKeys(const char *match = NULL);

  /** Write this Registry to an file.  The second parameter is an optional
   * header, each line of which must start with the '#" character  */
  void WriteToFile(const char *pathname, const char *header = NULL);

  /** Read this Registry from a file */
  void ReadFromFile(const char *pathname);

  /** Experimental */
  void SetFlagAddIfNotFound(bool yesno);

  /** Put an array into the registry */
  template <class T> void PutArray(unsigned int size,const T *array)
  {
    RemoveKeys("Element");
    Entry("ArraySize") << size;
    for(unsigned int i=0;i<size;i++)
      {
      Entry(Key("Element[%d]",i)) << array[i];
      }
  }

  /** Put an array into the registry */
  template <class T> void PutArray(const std::vector<T> &array)
  {
    RemoveKeys("Element");
    Entry("ArraySize") << array.size();
    for(unsigned int i=0;i<array.size();i++)
      {
      Entry(Key("Element[%d]",i)) << array[i];
      }
  }

  /** Get the size of an array */
  unsigned int GetArraySize() 
    { return Entry("ArraySize")[(unsigned int) 0]; }

  /** Get an array from the registry and store it in the vector passed in */
  template <class T> void GetArray(T *array, const T &defaultElement)
    {
    // Try reading the element count
    unsigned int size = Entry("ArraySize")[(unsigned int) 0];

    // Read element
    for(unsigned int i=0;i < size;i++)
      array[i] = Entry(Key("Element[%d]",i))[defaultElement];
    }

  /** Get an array from the registry */
  template <class T> std::vector<T> GetArray(const T &defaultElement)
  {
    // Try reading the element count
    unsigned int size = Entry("ArraySize")[(unsigned int) 0];

    // Initialize the result array
    std::vector<T> result(size,defaultElement);

    // Read element
    for(unsigned int i=0;i < size;i++)
      {
      result[i] = Entry(Key("Element[%d]",i))[defaultElement];
      }

    // Return the array
    return result;
  }

  /** An IO exception objects thrown by this class when reading from file*/
  class IOException : public StringType {
  public:
    IOException(const char *text) : StringType(text) {}
  };

  /** A syntax error exception */
  class SyntaxException : public StringType {
  public:
    SyntaxException(const char *text) : StringType(text) {}
  };

  Registry &operator =(const Registry &rsrc);

private:
  // Hashtable type definition
  typedef std::map<StringType, Registry *> FolderMapType;
  typedef std::map<StringType, RegistryValue> EntryMapType;

  // Commonly used hashtable iterators
  typedef FolderMapType::const_iterator FolderIterator;
  typedef EntryMapType::iterator EntryIterator;

  /** A hash table for the subfolders */
  FolderMapType m_FolderMap;

  /** A hash table for the registry values */
  EntryMapType m_EntryMap;

  /**
   * A flag as to whether keys and folders that are read and not found
   * should be created and populated with default values.
   */
  bool m_AddIfNotFound;

  /** Write this folder recursively to a stream */
  void Write(std::ostream &sout,const StringType &keyPrefix);

  /** Read this folder recursively from a stream, recording syntax errors */
  void Read(std::istream &sin, std::ostream &serr);

  /** Encode a string for writing to file */
  static StringType Encode(const StringType &input);

  /** Decode a string for writing to file */
  static StringType Decode(const StringType &input);
};

#endif
