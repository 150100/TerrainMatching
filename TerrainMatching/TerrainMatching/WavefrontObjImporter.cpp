///////////////////////////////////////////////////////////////////////////////
//
//  Copyright Leonardo Fischer 2011 - http://coderender.blogspot.com
//
//  Distributed under the licence available in the accompanying file 
//  LICENCE.txt. Please read it before use this code.
//
///////////////////////////////////////////////////////////////////////////////

#include "WavefrontObjImporter.h"
#include <fstream>
#include <sstream>

void WavefrontObjLoader::load( const std::string& objFile )
{
    vertices.clear();
    verticeCount=0;
    faces.clear();
    faceCount=0;

    std::ifstream stream;
    try
    {
        stream.open( objFile.c_str(), std::ios::in );
        if( stream )
        {
            // load the obj
            while( stream.good() )
            {
                char buffer[512];
                stream.getline(buffer, sizeof(buffer));
                std::stringstream str(buffer);
                std::string name;
                str >> name;
                if( name=="v" )
                {
                    Vector3f position;
                    str >> position.x >> position.y >> position.z;
                    vertices.push_back( position );
                    verticeCount++;
                }

                /* modified by alsrbbk */

                else if( name=="f" )
                {
                    unsigned int v1, v2, v3; // vertex geometry
                    unsigned int vt1, vt2, vt3; // vertex texture coordinate
                    unsigned int vn1, vn2, vn3; // vertex normal
                    v1=v2=v3=vt1=vt2=vt3=vn1=vn2=vn3=0;

                    char buff;
                    str >> v1;
                    str.get(buff);
                    if (buff == '/') {
                        str.get(buff);
                        if (buff == '/') {
                            // format of "f [v1]//[vn1] [v2]//[vn2] [v3]//[vn3]
                            str >> vn1 >> v2 >> buff >> buff >> vn2 >> v3 >> buff >> buff >> vn3;
                        }
                        else if (buff >= '0' && buff <= '9') {
                            str.unget();
                            str >> vt1;
                            str.get(buff);
                            if (buff == '/') {
                                // format of "f [v1]/[vt1]/[vn1] [v2]/[vt2]/[vn2] [v3]/[vt3]/[vn3]
                                str >> vn1 >> v2 >> buff >> vt2 >> buff >> vn2 >> v3 >> buff >> vt3 >> buff >> vn3;
                            }
                            else if (buff == ' ') {
                                // format of "f [v1]/[vt1] [v2]/[vt2] [v3]/[vt3]
                                str >> v2 >> buff >> vt2 >> v3 >> buff >> vt3;
                            }
                        }
                        else
                            throw cpp::Exception(".obj format does not match.");
                    }
                    else if(buff == ' ') {
                        // format of "f [v1] [v2] [v3]"
                        str >> v2 >> v3;
                    }
                    else
                        throw cpp::Exception(".obj format does not match.");

                    faces.push_back(v1-1);
                    faces.push_back(v2-1);
                    faces.push_back(v3-1);
                    faceCount++;
                }

                /* modified by alsrbbk end */
            }
            stream.close();
        }
        else
        {
            throw cpp::Exception("Can't open the file '"+objFile+"'");
        }
    }
    catch( const std::exception& e )
    {
        if( stream.is_open() )
            stream.close();
        std::cerr << e.what();
        throw;
    }
}
