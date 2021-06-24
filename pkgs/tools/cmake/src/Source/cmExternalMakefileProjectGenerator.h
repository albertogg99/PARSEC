/*=========================================================================

  Program:   CMake - Cross-Platform Makefile Generator
  Module:    $RCSfile: cmExternalMakefileProjectGenerator.h,v $
  Language:  C++
  Date:      $Date: 2012/03/29 17:21:08 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) 2007 Kitware, Inc., Insight Consortium.  All rights reserved.
  See Copyright.txt or http://www.cmake.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef cmExternalMakefileProjectGenerator_h
#define cmExternalMakefileProjectGenerator_h

#include "cmStandardIncludes.h"

#include "cmDocumentation.h"

class cmGlobalGenerator;

/** \class cmExternalMakefileProjectGenerator
 * \brief Base class for generators for "External Makefile based IDE projects".
 *
 * cmExternalMakefileProjectGenerator is a base class for generators
 * for "external makefile based projects", i.e. IDE projects which work 
 * an already existing makefiles.
 * See cmGlobalKdevelopGenerator as an example.
 * After the makefiles have been generated by one of the Makefile 
 * generators, the Generate() method is called and this generator
 * can iterate over the local generators and/or projects to produce the 
 * project files for the IDE.
 */
class cmExternalMakefileProjectGenerator
{
public:

  virtual ~cmExternalMakefileProjectGenerator() {}

  ///! Get the name for this generator.
  virtual const char* GetName() const = 0;
  /** Get the documentation entry for this generator.  */
  virtual void GetDocumentation(cmDocumentationEntry& entry, 
                                const char* fullName) const = 0;

  ///! set the global generator which will generate the makefiles
  virtual void SetGlobalGenerator(cmGlobalGenerator* generator)
                                           {this->GlobalGenerator = generator;}

  ///! Return the list of global generators supported by this extra generator
  const std::vector<std::string>& GetSupportedGlobalGenerators() const 
                                      {return this->SupportedGlobalGenerators;}

  ///! Get the name of the global generator for the given full name
  const char* GetGlobalGeneratorName(const char* fullName);
  /** Create a full name from the given global generator name and the
   * extra generator name
   */
  static std::string CreateFullGeneratorName(const char* globalGenerator, 
                                             const char* extraGenerator);

  ///! Generate the project files, the Makefiles have already been generated
  virtual void Generate() = 0;
protected:
  ///! Contains the names of the global generators support by this generator.
  std::vector<std::string> SupportedGlobalGenerators;
  ///! the global generator which creates the makefiles
  const cmGlobalGenerator* GlobalGenerator;
};

#endif