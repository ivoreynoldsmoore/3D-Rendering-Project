#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <string>
#include <iostream>
#include <ModelTriangle.h>
#include <TextureMap.h>
#include <map>
#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#include <RayTriangleIntersection.h>
//#define WIDTH 320
#define WIDTH 640
//#define HEIGHT 240
#define HEIGHT 480
#define SCALEFACTOR 0.17
#define MAXDEPTH 2

struct ModelObject {
    std::string name;
    std::vector<int> faces{};
    bool reflective;
    bool transparent;
    float refractiveIndex;

    ModelObject();
};

struct OBJSpecs {
    std::string filename;
    bool clockwise;
    float scaleFactor;
    glm::vec3 offset;
};

ModelObject::ModelObject() = default;

float magnitude(glm::vec3 a) {
    return sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2));
}

glm::vec3 normalise(glm::vec3 a) {
    float Mag = magnitude(a);
    glm::vec3 b(a[0]/Mag, a[1]/Mag, a[2]/Mag);
    return b;
}

float findArea(glm::vec3 pointA, glm::vec3 pointB, glm::vec3 pointC) {
    float a = glm::length((pointA-pointB));
    float b = glm::length((pointA-pointC));
    float c = glm::length((pointB-pointC));
    float s = (a+b+c)/2;
    float area = sqrt(glm::max((s*(s-a)*(s-b)*(s-c)),0.0f));
    return area;
}

float findArea(glm::vec2 pointA, glm::vec2 pointB, glm::vec2 pointC) {
    float a = glm::length((pointA-pointB));
    float b = glm::length((pointA-pointC));
    float c = glm::length((pointB-pointC));
    float s = (a+b+c)/2;
    float area = sqrt(glm::max((s*(s-a)*(s-b)*(s-c)),0.0f));
    return area;
}

float findArea(ModelTriangle triangle) {
    float a = magnitude((triangle.vertices[0]-triangle.vertices[1]));
    float b = magnitude((triangle.vertices[0]-triangle.vertices[2]));
    float c = magnitude((triangle.vertices[1]-triangle.vertices[2]));
    float s = (a+b+c)/2;
    float area = sqrt(glm::max((s*(s-a)*(s-b)*(s-c)),0.0f));
    return area;
}

TexturePoint findTexturePoint(ModelTriangle triangle, glm::vec3 point, bool debug) {
    glm::vec3 baryCoord;
    TexturePoint texturePoint;
    float triangleArea = findArea(triangle);
    baryCoord[0] = findArea(triangle.vertices[1], triangle.vertices[2], point)/triangleArea;
    baryCoord[1] = findArea(triangle.vertices[0], triangle.vertices[2], point)/triangleArea;
    baryCoord[2] = findArea(triangle.vertices[0], triangle.vertices[1], point)/triangleArea;
//    std::cout << "Bary Coord: (" << baryCoord[0] << ", " << baryCoord[1] << ", " << baryCoord[2] << ")" << '\n';
    texturePoint.x = round(baryCoord[0]*triangle.texturePoints[0].x+baryCoord[1]*triangle.texturePoints[1].x+baryCoord[2]*triangle.texturePoints[2].x);
    texturePoint.y = round(baryCoord[0]*triangle.texturePoints[0].y+baryCoord[1]*triangle.texturePoints[1].y+baryCoord[2]*triangle.texturePoints[2].y);
    if (debug) {
    }
    return texturePoint;
}

TexturePoint findTexturePoint(CanvasTriangle triangle, glm::vec2 point) {
    glm::vec3 baryCoord;
    TexturePoint texturePoint;
    glm::vec2 v0(triangle.vertices[0].x, triangle.vertices[0].y);
    glm::vec2 v1(triangle.vertices[1].x, triangle.vertices[1].y);
    glm::vec2 v2(triangle.vertices[2].x, triangle.vertices[2].y);

    baryCoord[0] = findArea(v1, v2, point);
    baryCoord[1] = findArea(v0, v2, point);
    baryCoord[2] = findArea(v0, v1, point);
    float triangleArea = baryCoord[0]+baryCoord[1]+baryCoord[2];
    baryCoord[0] = (baryCoord[0]/triangleArea);
    baryCoord[1] = (baryCoord[1]/triangleArea);
    baryCoord[2] = (baryCoord[2]/triangleArea);
    texturePoint.x = round(baryCoord[0]*triangle.texturePoints[0].x+baryCoord[1]*triangle.texturePoints[1].x+baryCoord[2]*triangle.texturePoints[2].x);
    texturePoint.y = round(baryCoord[0]*triangle.texturePoints[0].y+baryCoord[1]*triangle.texturePoints[1].y+baryCoord[2]*triangle.texturePoints[2].y);
    return texturePoint;
}

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
    float step = (to - from)/(numberOfValues-1);
    std::vector<float> floats;
    for (int i = 0; i < numberOfValues; i++){
        floats.push_back(from+(step*i));
    }
    return floats;
}

std::vector<glm::vec2> interpolate2DCoords(CanvasPoint from, CanvasPoint to, int numberOfValues) {
    std::vector<glm::vec2> vectors;
    std::vector<float> vals1 = interpolateSingleFloats(from.x, to.x, numberOfValues);
    std::vector<float> vals2 = interpolateSingleFloats(from.y, to.y, numberOfValues);
    for(int i = 0; i < numberOfValues; i++){
        vectors.push_back(glm::vec2(vals1[i], vals2[i]));
    }
    return vectors;
}

std::vector<glm::vec2> interpolate2DCoords(glm::vec2 from, glm::vec2 to, int numberOfValues) {
    std::vector<glm::vec2> vectors;
    std::vector<float> vals1 = interpolateSingleFloats(from[0], to[0], numberOfValues);
    std::vector<float> vals2 = interpolateSingleFloats(from[1], to[1], numberOfValues);
    for(int i = 0; i < numberOfValues; i++){
        vectors.push_back(glm::vec2(vals1[i], vals2[i]));
    }
    return vectors;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues){
    std::vector<glm::vec3> vectors;
    std::vector<float> vals1 = interpolateSingleFloats(from[0], to[0], numberOfValues);
    std::vector<float> vals2 = interpolateSingleFloats(from[1], to[1], numberOfValues);
    std::vector<float> vals3 = interpolateSingleFloats(from[2], to[2], numberOfValues);
    for(int i = 0; i < numberOfValues; i++){
        vectors.push_back(glm::vec3(vals1[i], vals2[i], vals3[i]));
    }
    return vectors;
}

std::vector<CanvasPoint> interpolateThreeElementValues(CanvasPoint from, CanvasPoint to, int numberOfValues) {
    std::vector<CanvasPoint> vectors;
    std::vector<float> vals1 = interpolateSingleFloats(from.x, to.x, numberOfValues);
    std::vector<float> vals2 = interpolateSingleFloats(from.y, to.y, numberOfValues);
    std::vector<float> vals3 = interpolateSingleFloats(from.depth, to.depth, numberOfValues);
    for(int i = 0; i < numberOfValues; i++){
        vectors.push_back(CanvasPoint(vals1[i], vals2[i], vals3[i]));
    }
    return vectors;
}

void importMTL(std::string filename, std::map<std::string, Colour>& colours, std::map<std::string, TextureMap>& textures) {
    std::ifstream file (filename);
    std::string values[3] = {};
    std::string materialStr;
    std::string ppmFilename;
    std::string label;

    while(!file.eof()) {
        file >> label;
//        std::cout << label << '\n';
        if (label.compare("newmtl") == 0) {
            file >> materialStr;
//            std::cout << "newmtl " << materialStr << '\n';
        }
        else if (label.compare("Kd") == 0) {
            file >> values[0] >> values[1] >> values[2];
            Colour colour(std::stof(values[0])*255.0,std::stof(values[1])*255.0,std::stof(values[2])*255.0);
            colour.alpha = 255.0;
            colours[materialStr] = colour;
//            std::cout << "Kd " << colours[materialStr].red << " " << colours[materialStr].green << " " << colours[materialStr].blue << '\n';
        }
        else if (label.compare("map_Kd") == 0) {
            file >> ppmFilename;
            TextureMap texture(ppmFilename);
            textures[materialStr] = texture;
//            std::cout << "map_Kd " << ppmFilename << '\n';
        }
    }
    file.close();
}

void importOBJ(OBJSpecs fileSpecs, std::vector<ModelTriangle>& faces, std::vector<ModelObject>& objects, std::map<std::string, TextureMap>& textures) {
    std::string values[3] = {};
    int textureIndexes[3] = {};
    std::string materialStr;
    std::string label;
    std::string mtlFilename;
    bool hasNormals = false;
    bool textured = false;

    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> vertexNormals;
    std::vector<glm::vec2> textureVertices;
    std::string obj;
    std::map<std::string, Colour> colours;

    std::ifstream file (fileSpecs.filename);
    file >> label >> mtlFilename;
    importMTL(mtlFilename, colours, textures);

//    std::cout << colours[materialStr].red << ", " << colours[materialStr].green << ", " << colours[materialStr].blue<< '\n';
    while(!file.eof()) {
        file >> label;
        if (label.compare("usemtl") == 0) {
            file >> materialStr;
//            std::cout << "usemtl " << materialStr << '\n';
        }
        else if (label.compare("o") == 0) {
            file >> obj;
            ModelObject temp;
            temp.name = obj;
            temp.reflective = false;
            temp.transparent = false;
            textured = false;
            hasNormals = false;
            temp.refractiveIndex = 1;
            objects.push_back(temp);
//            std::cout << "o " << obj << '\n';
        }
        else if (label.compare("rl") == 0) {
            objects[objects.size()-1].reflective = true;
        }
        else if (label.compare("tp") == 0) {
            objects[objects.size()-1].transparent = true;
        }
        else if (label.compare("ri") == 0) {
            file >> values[0];
            objects[objects.size()-1].refractiveIndex = std::stof(values[0]);
        }
        else if (label.compare( "v") == 0) {
            file >> values[0] >> values[1] >> values[2];
//            std::cout << values[0] << '\t' << values[1] << '\t'  << values[2] << '\n';
            glm::vec3 temp = glm::vec3(std::stof(values[0])*fileSpecs.scaleFactor, std::stof(values[1])*fileSpecs.scaleFactor, std::stof(values[2])*fileSpecs.scaleFactor)+fileSpecs.offset;
            vertices.push_back(temp);
//            std::cout << "v " << values[0] << " " << values[1] << " " << values[2] << '\n';
        }
        else if (label.compare( "vt") == 0) {
            file >> values[0] >> values[1];
            textured = true;
//            std::cout << values[0] << '\t' << values[1] << '\n';
            glm::vec2 temp = glm::vec2(std::stof(values[0])*(textures[materialStr].width-1), std::stof(values[1])*(textures[materialStr].height-1));
            textureVertices.push_back(temp);
//            std::cout << "vt " << values[0] << " " << values[1] << " " << values[2] << '\n';
        }
        else if (label.compare( "vn") == 0) {
            hasNormals = true;
            file >> values[0] >> values[1] >> values[2];
            glm::vec3 temp = glm::normalize(glm::vec3(std::stof(values[0]), std::stof(values[1]), std::stof(values[2])));
            vertexNormals.push_back(temp);
//            std::cout << "vn " << values[0] << " " << values[1] << " " << values[2] << '\n';
        }
        else if (label.compare("f") == 0) {
            file >> values[0] >> values[1] >> values[2];
            std::string sep = "/";
            ModelTriangle temp;
            temp.textured = textured;
            for (int i = 0; i < 3; i++) {
                if (temp.textured) {
                    textureIndexes[i] = std::stoi(values[i].substr(values[i].find(sep)+1));
                    temp.texturePoints[i] = TexturePoint(textureVertices[textureIndexes[i]-1][0], textureVertices[textureIndexes[i]-1][1]);
                }
                values[i] = values[i].substr(0, values[i].find(sep));
                temp.vertices[i] = vertices[std::stoi(values[i])-1];
                if (hasNormals){
                    temp.vertexNormals[i] = vertexNormals[std::stoi(values[i])-1];
                }
            }
            if (fileSpecs.clockwise) {
                temp.normal = glm::normalize(glm::cross(temp.vertices[1]-temp.vertices[0], temp.vertices[2]-temp.vertices[0]));
            }
            else {
                temp.normal = glm::normalize(glm::cross(temp.vertices[2]-temp.vertices[0], temp.vertices[1]-temp.vertices[0]));
            }
            temp.objectIndex = (objects.size()-1);
            temp.colour = colours[materialStr];
            temp.colour.name = materialStr;
            temp.hasVertexNormals = hasNormals;
            if (temp.textured) {
//                std::cout << "Texture: " << temp.colour.name << '\n';
//                std::cout << textureIndexes[0] << ", " << textureIndexes[1] << ", " << textureIndexes[2] << '\n';
//                std::cout << "("<< temp.texturePoints[0] << "), (" << temp.texturePoints[1] << "), (" << temp.texturePoints[2] << ") \n";
//                std::cout << faces.size()+1 << '\n';
            }
            faces.push_back(temp);
//            std::cout << "f " << values[0] << " " << values[1] << " " << values[2] << '\n';
        }
    }
    file.close();
}

void drawline(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour){
    float dx = from.x - to.x;
    float dy = from.y - to.y;
    int numOfVals = std::round(std::max(fabs(dx), fabs(dy)));
    uint32_t colInt = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
    std::vector<CanvasPoint> line = interpolateThreeElementValues(from, to, numOfVals);
    for (int i = 0; i < line.size(); i++){
        if (std::round(line[i].x) > 0 && std::round(line[i].x) < WIDTH && std::round(line[i].y) > 0 && std::round(line[i].y) < HEIGHT) {
            window.setPixelColour(std::round(line[i].x), std::round(line[i].y), colInt);
//            std::cout << "Set Pixel" << '\n';
        }
    }
}

void drawTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour){
    drawline(window, triangle.v0(), triangle.v1(), colour);
    drawline(window, triangle.v1(),triangle.v2(), colour);
    drawline(window, triangle.v2(), triangle.v0(), colour);
}

void fillTriangleHalf(DrawingWindow &window, CanvasTriangle triangle, std::map< std::string, std::vector< std::vector< Colour > > >& textures, std::vector<CanvasPoint> sideA, std::vector<CanvasPoint> sideB, std::vector<std::vector<float>>& depthBuffer) {
    Colour colour = triangle.colour;
    TexturePoint texturePoint;
    int start = 0;
    int end = 0;
    for (int y = std::round(sideA[0].y); y < sideA[sideA.size()-1].y; y++) {
//        std::cout << '\n' << "Current y: " << y << '\n';
        for (int i = 0; i < int(sideA.size()); i++) {
            if (std::round(sideA[i].y) == y) {
                start = i;
                break;
            }
        }
        for (int i = 0; i < sideB.size(); i++) {
            if (std::round(sideB[i].y) == y) {
                end = i;
                break;
            }
        }

        std::vector<float> depthRow = interpolateSingleFloats(sideA[start].depth, sideB[end].depth, fabs(std::round(sideA[start].x)-std::round(sideB[end].x))+1);
        std::vector<float> row = interpolateSingleFloats(std::round(sideA[start].x), std::round(sideB[end].x), fabs(std::round(sideA[start].x)-std::round(sideB[end].x))+1);
        for (int x = 0; x < row.size(); x++) {
            if (row[x] > 0 && row[x] < WIDTH && y > 0 && y < HEIGHT) {
//                std::cout << "Current x: " << row[x] << '\n';
                if ((depthBuffer[row[x]][y] == 0.0) || (depthBuffer[row[x]][y] > (1/depthRow[x]) && depthBuffer[row[x]][y] < 0.0)) {
                    if (triangle.textured) {
                        texturePoint = findTexturePoint(triangle, glm::vec2(row[x],y));
                        std::vector< std::vector< Colour > >& texture = textures[triangle.colour.name];
                        colour = texture[texturePoint.x][texturePoint.y];
//                        colour.alpha = texture[texturePoint.x][texturePoint.y].alpha;
//                        colour.red = texture[texturePoint.x][texturePoint.y].red;
//                        colour.blue = texture[texturePoint.x][texturePoint.y].blue;
//                        colour.green = texture[texturePoint.x][texturePoint.y].green;
//                        std::vector< std::vector< Colour > >& texture = textures[rayData.intersectedTriangle.colour.name];
                    }
                    uint32_t colInt = static_cast<int>(colour.alpha << 24) + static_cast<int>(colour.red << 16) + static_cast<int>(colour.green << 8) + static_cast<int>(colour.blue);
                    window.setPixelColour(row[x], std::round(y), colInt);
//                    std::cout << "Set Pixel" << '\n';
                    depthBuffer[row[x]][y] = 1/depthRow[x];
//                    std::cout << "Assigned z at: (" << row[x] << "," << y << ")" << '\n';
                }
            }
            else
            {
//                std::cout << "(" << row[x] << ", " << y << ")" <<  '\n';
            }
        }
    }
}

void fillTriangle(DrawingWindow &window, CanvasTriangle triangle, std::map< std::string, std::vector< std::vector< Colour > > >& textures, std::vector<std::vector<float>>& depthBuffer) {
    if (triangle.vertices[0].y > triangle.vertices[1].y){
        std::swap(triangle.vertices[0],triangle.vertices[1]);
        std::swap(triangle.texturePoints[0],triangle.texturePoints[1]);
    }
    if (triangle.vertices[1].y > triangle.vertices[2].y){
        std::swap(triangle.vertices[1],triangle.vertices[2]);
        std::swap(triangle.texturePoints[1],triangle.texturePoints[2]);
    }
    if(triangle.vertices[0].y > triangle.vertices[1].y){
        std::swap(triangle.vertices[0],triangle.vertices[1]);
        std::swap(triangle.texturePoints[0],triangle.texturePoints[1]);
}

//    std::cout << "Top point: " << triangle.vertices[0] << '\n';
//    std::cout << "Middle point: " << triangle.vertices[1] << '\n';
//    std::cout << "Bottom point: " << triangle.vertices[2] << '\n';

    int numOfVals = std::max(fabs(triangle.vertices[0].y-triangle.vertices[2].y), fabs(triangle.vertices[0].x-triangle.vertices[2].x));
    std::vector<CanvasPoint> sideT = interpolateThreeElementValues(triangle.vertices[0], triangle.vertices[2], numOfVals);

    std::vector<CanvasPoint> side1;
    std::vector<CanvasPoint> side4;
    for (int i = 0; i < sideT.size(); i++) {
//        std::cout << "Side T i: " << sideT[i][1] << '\n';
//        std::cout << "Current y: " << std::round(triangle.vertices[1].y) << '\n';
        if (sideT[i].y < std::round(triangle.vertices[1].y)) {
            side1.push_back(sideT[i]);
        }
        else if (sideT[i].y >= std::round(triangle.vertices[1].y)) {
            side4.push_back(sideT[i]);
        }
    }
    numOfVals = std::max(fabs(triangle.vertices[0].y-triangle.vertices[1].y), fabs(triangle.vertices[0].x-triangle.vertices[1].x));
    std::vector<CanvasPoint> side2;
    if (numOfVals != 1) {
        side2 = interpolateThreeElementValues(triangle.vertices[0], triangle.vertices[1], numOfVals);
    }
    else {
        side2 = {CanvasPoint((triangle.vertices[0].x+triangle.vertices[1].x),(triangle.vertices[0].y+triangle.vertices[1].y))};
    }

    numOfVals = std::max(fabs(triangle.vertices[1].y-triangle.vertices[2].y), fabs(triangle.vertices[1].x-triangle.vertices[2].x));
    std::vector<CanvasPoint> side3;
    if (numOfVals != 1) {
        side3 = interpolateThreeElementValues(triangle.vertices[1], triangle.vertices[2], numOfVals);
    }
    else {
        side3 = {CanvasPoint((triangle.vertices[1].x+triangle.vertices[2].x),(triangle.vertices[1].y+triangle.vertices[2].y))};
    }
    if (side1.size() > 0) {
//        std::cout << "First Triangle" << '\n';
        CanvasTriangle triangle1;
        triangle1.vertices = {triangle.vertices[0], triangle.vertices[1], side1.back()};
        triangle1.texturePoints = {triangle.texturePoints[0], triangle.texturePoints[1]};
        fillTriangleHalf(window, triangle, textures, side1, side2, depthBuffer);
    }

    if (side4.size() > 0) {
//        std::cout << "Second Triangle" << '\n';
        fillTriangleHalf(window, triangle, textures, side3, side4, depthBuffer);
    }
//    std::cout << "Triangles Done" << '\n';
}

bool getClosestIntersection(std::vector<ModelTriangle>& faces, std::vector<ModelObject>& objects, glm::vec3 rayOrigin, glm::vec3 rayDirection, RayTriangleIntersection& rayData, bool light, bool debug = false) {
    glm::vec3 e0;
    glm::vec3 e1;
    glm::vec3 SPVector;
    glm::vec3 possibleSolution;
    glm::mat3 DEMatrix;
    glm::vec3 closestTri(-FLT_MAX, 0, 0);
    std::size_t closestTriIndex = 0;
    bool found = false;
    for (int i = 0; i < faces.size(); i++){
        RayTriangleIntersection result;
        e0 = faces[i].vertices[1] - faces[i].vertices[0];
        e1 = faces[i].vertices[2] - faces[i].vertices[0];
        SPVector = rayOrigin - faces[i].vertices[0];
        DEMatrix = glm::mat3(-rayDirection, e0, e1);
        possibleSolution = glm::inverse(DEMatrix) * SPVector;
        if (debug){
//            std::cout << "Face Index: " << i << '\n';
//            std::cout << "Possible Solution: (" << possibleSolution[0] << ", " << possibleSolution[1] << ", " << possibleSolution[2] << ")" <<  '\n';
        }
        if ((possibleSolution[0] > closestTri[0])
        && possibleSolution[0] < 0
        && possibleSolution[1] <= 1 && possibleSolution[1] >= 0
        && possibleSolution[2] <= 1 && possibleSolution[2] >= 0
        && (possibleSolution[1] + possibleSolution[2]) <= 1
        && !(light && objects[faces[i].objectIndex].transparent)) {
            closestTri = possibleSolution;
            closestTriIndex = i;
//            std::cout << faces[i] << '\n';
            found = true;
        }
    }
    if (found) {
        if (debug) {
//            std::cout << "Accepted Solution: (" << closestTri[0] << ", " << closestTri[1] << ", " << closestTri[2] << ")" <<  '\n';
//            std::cout << "Chosen Face Index: " << closestTriIndex << '\n';
        }
        rayData = RayTriangleIntersection((rayOrigin + closestTri[0]*rayDirection), fabs(closestTri[0]), faces[closestTriIndex], closestTriIndex);
    }
    return found;
}

glm::vec3 modCoordToCamCoord(glm::vec3 cameraPosition, glm::vec3 vertexPosition) {
//    std::cout << "Coordinate Conversion: (" << vertexPosition[0] << ", " << vertexPosition[1] << ", " << vertexPosition[2] << ") \n";
    glm::vec3 newVertexPosition;
    newVertexPosition = (vertexPosition - cameraPosition);
//    std::cout << "New Coordinate: (" << newVertexPosition[0] << ", " << newVertexPosition[1] << ", " << newVertexPosition[2] << ") \n";
    return newVertexPosition;
}

glm::vec2 getCanvasIntersectionPoint(glm::vec3 vertexPosition, float focalLength) {
    float newX = (WIDTH/2)-((focalLength/vertexPosition[2])*vertexPosition[0]*256);
    float newY = (HEIGHT/2)+((focalLength/vertexPosition[2])*vertexPosition[1]*256);
    glm::vec2 newVertexPosition(newX, newY);
    return newVertexPosition;
}

void lookAt(glm::vec3 point, glm::vec3 cameraPosition, glm::mat3& cameraOrientation) {
    glm::vec3 forward = glm::normalize(-modCoordToCamCoord(cameraPosition, point));
    glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0,1,0), cameraOrientation[2]));
    glm::vec3 up = glm::normalize(glm::cross(cameraOrientation[2], cameraOrientation[0]));
    cameraOrientation[2] = forward;
    cameraOrientation[0] = right;
    cameraOrientation[1] = up;
}

std::vector<CanvasTriangle> generateImagePlane(std::vector<ModelTriangle>& faces, float focalLength, glm::vec3 cameraPosition, glm::mat3 cameraOrientation) {
    std::vector<CanvasTriangle> imagePlane = {};
    CanvasTriangle newTriangle;
    glm::vec3 camRelativeVertex1;
    glm::vec3 camRelativeVertex2;
    glm::vec3 camRelativeVertex3;

    glm::vec2 newV0;
    glm::vec2 newV1;
    glm::vec2 newV2;

    for (int i = 0; i < faces.size(); i++) {

        camRelativeVertex1 = modCoordToCamCoord(cameraPosition, faces[i].vertices[0])*cameraOrientation;
        camRelativeVertex2 = modCoordToCamCoord(cameraPosition, faces[i].vertices[1])*cameraOrientation;
        camRelativeVertex3 = modCoordToCamCoord(cameraPosition, faces[i].vertices[2])*cameraOrientation;


        newV0 = getCanvasIntersectionPoint(camRelativeVertex1, focalLength);
        newV1 = getCanvasIntersectionPoint(camRelativeVertex2, focalLength);
        newV2 = getCanvasIntersectionPoint(camRelativeVertex3, focalLength);

        newTriangle = CanvasTriangle();
        newTriangle.vertices[0] = CanvasPoint(newV0[0], newV0[1], camRelativeVertex1[2]);
        newTriangle.vertices[1] = CanvasPoint(newV1[0], newV1[1], camRelativeVertex2[2]);
        newTriangle.vertices[2] = CanvasPoint(newV2[0], newV2[1], camRelativeVertex3[2]);
        newTriangle.textured = faces[i].textured;
        newTriangle.texturePoints = faces[i].texturePoints;
        newTriangle.colour = faces[i].colour;
        imagePlane.push_back(newTriangle);
    }
    return imagePlane;
}

void drawRasterisedScene(DrawingWindow &window, std::vector<ModelTriangle>& faces, std::map< std::string, std::vector< std::vector< Colour > > >& textures, std::vector<std::vector<float>> depthBuffer, float focalLength, glm::vec3& cameraPosition, glm::mat3& cameraOrientation, bool orbit) {
    if (orbit) {
        float rotationAngle = 0.01;
        cameraPosition = cameraPosition*glm::mat3(
                cos(rotationAngle), 0, sin(rotationAngle),
                0, 1, 0,
                -sin(rotationAngle), 0,  cos(rotationAngle)
        );
        lookAt(glm::vec3(0,0,0), cameraPosition, cameraOrientation);
    }
    window.clearPixels();
    std::vector<CanvasTriangle> imagePlane = generateImagePlane(faces, focalLength, cameraPosition, cameraOrientation);
    for (int i = 0; i < imagePlane.size(); i++) {
        fillTriangle(window, imagePlane[i], textures, depthBuffer);
    }
}

void drawWireframe(DrawingWindow &window, std::vector<ModelTriangle>& faces, float focalLength, glm::vec3& cameraPosition, glm::mat3& cameraOrientation, bool orbit) {
    if (orbit) {
        float rotationAngle = 0.01;
        cameraPosition = cameraPosition*glm::mat3(
                cos(rotationAngle), 0, sin(rotationAngle),
                0, 1, 0,
                -sin(rotationAngle), 0,  cos(rotationAngle)
        );
        lookAt(glm::vec3(0,0,0), cameraPosition, cameraOrientation);
    }
    window.clearPixels();
    std::vector<CanvasTriangle> imagePlane = generateImagePlane(faces, focalLength, cameraPosition, cameraOrientation);
    for (int i = 0; i < imagePlane.size(); i++) {
        drawTriangle(window, imagePlane[i], Colour(255,255,255));
    }
}

glm::vec3 getDirectionVector(glm::vec2 vertexPosition, float focalLength) {
    float newX = (((WIDTH/2)-vertexPosition[0])/(256*focalLength));
    float newY = ((vertexPosition[1]-(HEIGHT/2))/(256*focalLength));
    float newZ = 1;
    glm::vec3 newVertexPosition(newX, newY, newZ);
    return glm::normalize(newVertexPosition);
}

glm::vec3 findNormal(ModelTriangle triangle, glm::vec3 point) {
    glm::vec3 baryCoord;

    float triangleArea = findArea(triangle);
    baryCoord[0] = findArea(triangle.vertices[1], triangle.vertices[2], point)/triangleArea;
    baryCoord[1] = findArea(triangle.vertices[0], triangle.vertices[2], point)/triangleArea;
    baryCoord[2] = findArea(triangle.vertices[0], triangle.vertices[1], point)/triangleArea;
//    std::cout << "Bary Coord: (" << baryCoord[0] << ", " << baryCoord[1] << ", " << baryCoord[2] << ")" << '\n';
    glm::vec3 normal = glm::normalize(baryCoord[0]*triangle.vertexNormals[0]+baryCoord[1]*triangle.vertexNormals[1]+baryCoord[2]*triangle.vertexNormals[2]);
    return normal;
}

glm::vec3 reflectRay(glm::vec3 incidentRay, glm::vec3 normal) {
    glm::vec3 reflectedRay = glm::normalize(incidentRay-2.0f*normal*glm::dot(incidentRay, normal));
    return reflectedRay;
}

glm::vec3 refractRay(glm::vec3 incidentRay, glm::vec3 normal, float n1, float n2, bool debug = false) {
    float c =  glm::dot(-normal, incidentRay);
    if (c < 0){
        float temp = n1;
        n1 = n2;
        n2 = temp;
        c = fabs(c);
    }
    float r = n1/n2;
    glm::vec3 refractedRay = r*incidentRay+float(r*c-sqrt(glm::max(float(1.0f-pow(r,2)*(1.0f-pow(c,2))), 0.0f)))*normal;
    if (debug) {
//        std::cout << "c: " << c << '\n';
//        std::cout << "r: " << r << '\n';
    }
    return refractedRay;
}

float calculateReflectivity(glm::vec3 incidentRay, glm::vec3 normal, float n1, float n2, bool debug = false) {
    float reflectivity;
    float theta =  glm::dot(-normal, incidentRay);
    if (theta < 0){
        float temp = n1;
        n1 = n2;
        n2 = temp;
        theta = fabs(theta);
    }
    float r0 = pow(((n1-n2)/(n1+n2)),2);
    reflectivity = r0 + (1-r0)*pow((1-theta),5);
    if (debug) {
//        std::cout << "r0: " << r0 << '\n';
//        std::cout << "Angle: " << theta << '\n';
    }
    return reflectivity;
}

float applyLighting(glm::vec3 point, float distance, glm::vec3 incidentRay, glm::vec3 normal,  glm::vec3 cameraPosition, bool debug = false){
    float brightness = 1;

    float proxBrightness = 2/(3.14*4*pow(distance/4, 2));
    proxBrightness = glm::clamp(proxBrightness,0.0f,1.0f);

    float aoiBrightness = glm::dot(normal, -incidentRay)*1.0f;
    aoiBrightness = glm::clamp(aoiBrightness,0.0f,1.0f);

    glm::vec3 reflectedRay = reflectRay(-incidentRay, normal);
    glm::vec3 view = -glm::normalize(cameraPosition - point);
    float specBrightness = pow(glm::max(glm::dot(glm::normalize(reflectedRay), view),0.0f),256)*5;
    specBrightness = glm::clamp(specBrightness,0.0f,1.0f);

    brightness = (proxBrightness*aoiBrightness)+specBrightness;
//    brightness = specBrightness;
    brightness = glm::clamp(brightness,0.0f,1.0f);
    if (debug){
//        std::cout << "Incident Ray: (" << -incidentRay[0] << ", " << -incidentRay[1] << ", " << -incidentRay[2] << ")" << '\n';
//        std::cout << "Proximity Lighting: " << proxBrightness <<'\n';
//        std::cout << "Normal: (" << normal[0] << ", " << normal[1] << ", " << normal[2] << ")" << '\n';
//        std::cout << "AoI Lighting: " << aoiBrightness <<'\n';
//        std::cout << "Reflected Ray: (" << reflectedRay[0] << ", " << reflectedRay[1] << ", " << reflectedRay[2] << ")" << '\n';
//        std::cout << "View: (" << view[0] << ", " << view[1] << ", " << view[2] << ")" << '\n';
//        std::cout << "Specular Lighting: " << specBrightness <<'\n';
//        std::cout << "Final Brightness: " << brightness << '\n';
    }
    return brightness;
}

Colour evaluateRay(std::vector<ModelTriangle>& faces, std::vector<ModelObject>& objects, std::vector<glm::vec3>& lights,
                   std::map< std::string, std::vector< std::vector< Colour > > >& textures, float focalLength,
                   glm::vec3 cameraPosition, glm::vec3 origin, glm::vec3 directionVector, int depth, bool debug = false) {
    Colour colour;
    RayTriangleIntersection rayData;
    RayTriangleIntersection lightData;
    float reflectivity;
    glm::vec3 normal;
    glm::vec3 lightVector;
    TexturePoint texturePoint;
    bool lightHit = false;
    float brightness = 1.0f;
    float totalBrightness = 0.0f;

    float hit = getClosestIntersection(faces, objects, origin, directionVector, rayData, false, debug);
    if (hit) {
        if (rayData.intersectedTriangle.hasVertexNormals) {
            normal = findNormal(rayData.intersectedTriangle, rayData.intersectionPoint);
        }
        else {
            normal = rayData.intersectedTriangle.normal;
        }
        if (objects[rayData.intersectedTriangle.objectIndex].transparent) {
            glm::vec3 refractedVector = refractRay(directionVector, normal, 1, objects[rayData.intersectedTriangle.objectIndex].refractiveIndex, debug);
            if (glm::dot(normal, directionVector) < 0){
                normal = -normal;
            }
            if (debug) {
                std::cout << "Refract" << '\n';
//                std::cout << "Refracted Vector: (" << refractedVector[0] << ", " << refractedVector[1] << ", " << refractedVector[2] << ")" << '\n';
            }
            colour = evaluateRay(faces, objects, lights, textures, focalLength, cameraPosition, rayData.intersectionPoint-(0.01f*normal), refractedVector, depth, debug);
            if (debug) {
                std::cout << "Refract return" << '\n';
            }
        }
        else {
            for (int i = 0; i < lights.size(); i++) {
                brightness = 0.0f;
                lightVector = (rayData.intersectionPoint+(0.01f*normal))-lights[i];
                lightHit = getClosestIntersection(faces, objects, rayData.intersectionPoint+(0.01f*normal), glm::normalize(lightVector), lightData, true);
                lightData.distanceFromCamera = magnitude(lightData.intersectionPoint-rayData.intersectionPoint);
                if (debug){
//                    std::cout << "Light distance: " << lightData.distanceFromCamera << '\n';
                }
                lightData.distanceFromCamera = std::min(lightData.distanceFromCamera, magnitude(lightVector));
                if (!lightHit){
                    lightData.distanceFromCamera = magnitude(lightVector);
                }
                if (lightData.distanceFromCamera == magnitude(lightVector)) {
                    brightness = applyLighting(rayData.intersectionPoint, lightData.distanceFromCamera, glm::normalize(lightVector), normal, cameraPosition, debug);
                }
                if (debug) {
//                    std::cout << "Light Point " << i << ": (" << lights[i][0] << ", " << lights[i][1] << ", " << lights[i][2] << ")" << '\n';
//                    std::cout << "Brightness for light " << i << ": " << brightness << '\n';
                }
                totalBrightness += brightness/lights.size();
            }
            totalBrightness = std::max(totalBrightness, 0.2f);
            if (rayData.intersectedTriangle.textured) {
                texturePoint = findTexturePoint(rayData.intersectedTriangle, rayData.intersectionPoint, debug);
                std::vector< std::vector< Colour > >& texture = textures[rayData.intersectedTriangle.colour.name];
                colour.alpha = texture[texturePoint.x][texturePoint.y].alpha;
                colour.red = texture[texturePoint.x][texturePoint.y].red;
                colour.blue = texture[texturePoint.x][texturePoint.y].blue;
                colour.green = texture[texturePoint.x][texturePoint.y].green;
            }
            else {
                colour.red = rayData.intersectedTriangle.colour.red;
                colour.blue = rayData.intersectedTriangle.colour.blue;
                colour.green = rayData.intersectedTriangle.colour.green;
            }
            colour.red = std::min(colour.red*totalBrightness, 255.0f);
            colour.blue = std::min(colour.blue*totalBrightness, 255.0f);
            colour.green = std::min(colour.green*totalBrightness, 255.0f);
        }
        if (objects[rayData.intersectedTriangle.objectIndex].reflective && depth < MAXDEPTH) {
            glm::vec3 reflectedVector = reflectRay(directionVector, normal);
//            std::cout << "Ray Vector: (" << directionVector[0] << ", " << directionVector[1] << ", " << directionVector[2] << ")" << '\n';
            reflectivity = calculateReflectivity(glm::normalize(directionVector), normal, 1, objects[rayData.intersectedTriangle.objectIndex].refractiveIndex, debug);
//            reflectivity = 0;
            if (debug) {
                std::cout << "Reflect" << '\n';
//                std::cout << "Reflectivity: " << reflectivity << '\n';
            }
            Colour reflectedColour = evaluateRay(faces, objects, lights, textures, focalLength, cameraPosition, rayData.intersectionPoint+(0.01f*normal), reflectedVector, depth+1, debug);
            if (debug) {
                std::cout << "Reflect Return" << '\n';
            }
            colour.red = ((colour.red*(1-reflectivity)) + (reflectedColour.red*reflectivity));
            colour.blue = ((colour.blue*(1-reflectivity)) + (reflectedColour.blue*reflectivity));
            colour.green = ((colour.green*(1-reflectivity)) + (reflectedColour.green*reflectivity));
        }
    }
    if (debug){
        std::cout << "Intersection Point: (" << rayData.intersectionPoint[0] << ", " << rayData.intersectionPoint[1] << ", " << rayData.intersectionPoint[2] << ")" << '\n';
//        std::cout << "Light distance: " << lightData.distanceFromCamera << '\n';
//        std::cout << "Actual Light Vector: (" << lightVector[0] << ", " << lightVector[1] << ", " << lightVector[2] << ")" << '\n';
//        std::cout << "Normal: (" << normal[0] << ", " << normal[1] << ", " << normal[2] << ")" << '\n';
//        std::cout << "Camera Position: (" << cameraPosition[0] << ", " << cameraPosition[1] << ", " << cameraPosition[2] << ")" << '\n';
//        std::cout << "Face Index: " << rayData.triangleIndex << '\n';
//        std::cout << "Light Hit: " << lightHit << '\n';
//        std::cout << "Hit: " << hit << '\n';
//        std::cout << "Depth: " << depth << '\n';
//        std::cout << "Colour: " << colour << '\n';
//        std::cout << "Texture Point: (" << texturePoint.x << ", " << texturePoint.y << ") \n";
//        std::cout << "Ray Vector: (" << directionVector[0] << ", " << directionVector[1] << ", " << directionVector[2] << ")" << '\n';
//        std::cout << "Light Vector: (" << rayData.intersectionPoint[0]-lightData.intersectionPoint[0] << ", " << rayData.intersectionPoint[1]-lightData.intersectionPoint[1] << ", " << rayData.intersectionPoint[2]-lightData.intersectionPoint[2] << ")" << '\n';
//        std::cout << "Light Face Index: " << lightData.triangleIndex << '\n';
//        std::cout << "Light Intersection Point: (" << lightData.intersectionPoint[0] << ", " << lightData.intersectionPoint[1] << ", " << lightData.intersectionPoint[2] << ")" << '\n';
//        std::cout << "Actual Light distance: " << magnitude(lightVector) << '\n';
//        std::cout << "Brightness: " << totalBrightness << '\n';
//        std::cout << "Ray data: " << '\n' << rayData << '\n';
//        std::cout << '\n';
    }
    return colour;
}

void drawRayTracing(DrawingWindow &window, std::vector<ModelTriangle>& faces, std::vector<ModelObject>& objects, std::vector<glm::vec3> lights, std::map< std::string, std::vector< std::vector< Colour > > >& textures, float focalLength, glm::vec3& cameraPosition, glm::mat3& cameraOrientation, bool orbit) {
    if (orbit) {
        float rotationAngle = 0.1;
        cameraPosition = cameraPosition*glm::mat3(
                cos(rotationAngle), 0, sin(rotationAngle),
                0, 1, 0,
                -sin(rotationAngle), 0,  cos(rotationAngle)
        );
        lookAt(glm::vec3(0,0,0), cameraPosition, cameraOrientation);
        std::cout << "(" << cameraPosition[0] << ", " << cameraPosition[1] << ", " << cameraPosition[2] << ")" <<  '\n';
        std::cout << "(" << cameraOrientation[0][0] << ", " << cameraOrientation[1][0] << ", " << cameraOrientation[2][0] << ")" <<  '\n';
        std::cout << "(" << cameraOrientation[0][1] << ", " << cameraOrientation[1][1] << ", " << cameraOrientation[2][1] << ")" <<  '\n';
        std::cout << "(" << cameraOrientation[0][2] << ", " << cameraOrientation[1][2] << ", " << cameraOrientation[2][2] << ")" <<  '\n' <<  '\n';
    }
    window.clearPixels();
    uint32_t colInt;
    Colour colour;
    glm::vec3 directionVector;

    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
//            std::cout << "Coordinates: " << x << ", " << y << '\n';
            directionVector = getDirectionVector(glm::vec2(x, y), focalLength)*glm::transpose(cameraOrientation);
            colour = evaluateRay(faces, objects, lights, textures, focalLength, cameraPosition, cameraPosition, directionVector, 0);
//            if (colour.red != 0 && colour.green != 0 && colour.blue != 0) {
//                std::cout << "Final Colour: " << colour << '\n';
//            }
            colInt = (static_cast<int>(std::round(255)) << 24)
                     + (static_cast<int>(std::round(colour.red)) << 16)
                     + (static_cast<int>(std::round(colour.green)) << 8)
                     + static_cast<int>(std::round(colour.blue));
            window.setPixelColour(x, y, colInt);
        }
    }
    colInt = (255 << 24)
             + (255 << 16)
             + (255 << 8)
             + 255;
//    for (int i = 0; i < lights.size(); i++) {
//        glm::vec2 lightSpot = getCanvasIntersectionPoint(modCoordToCamCoord(cameraPosition, lights[i])*cameraOrientation, focalLength);
//        if (lightSpot[0] > 0 && lightSpot[0] < WIDTH && lightSpot[1] > 0 && lightSpot[1] < HEIGHT) {
//            window.setPixelColour(lightSpot[0], lightSpot[1], colInt);
//        }
//    }
}

void handleEvent(SDL_Event event, DrawingWindow &window, std::vector<ModelTriangle>& faces, std::vector<ModelObject>& objects,
                 std::vector<glm::vec3>& lights, std::map< std::string, std::vector< std::vector< Colour > > >& textures, float focalLength,
                 glm::vec3& cameraPosition, glm::mat3& cameraOrientation, bool& orbit, std::string& renderType) {
    float rotationAngle = 0.05;
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_w) {
            cameraPosition += glm::vec3(0, 0, -0.1)*cameraOrientation;
        }
        else if (event.key.keysym.sym == SDLK_s) {
            cameraPosition += glm::vec3(0, 0, 0.1)*cameraOrientation;
        }
        else if (event.key.keysym.sym == SDLK_a) {
            cameraPosition += glm::vec3(-0.1, 0, 0)*cameraOrientation;
        }
        else if (event.key.keysym.sym == SDLK_d) {
            cameraPosition += glm::vec3(0.1, 0, 0)*cameraOrientation;
        }
        else if (event.key.keysym.sym == SDLK_1 && renderType.compare("wireframe") != 0) {
            renderType  = "wireframe";
            std::cout << "Render Type: " << renderType << '\n';
        }
        else if (event.key.keysym.sym == SDLK_2 && renderType.compare("rasterise") != 0) {
            renderType  = "rasterise";
            std::cout << "Render Type: " << renderType << '\n';
        }
        else if (event.key.keysym.sym == SDLK_3 && renderType.compare("raytracing") != 0) {
            renderType  = "raytracing";
            std::cout << "Render Type: " << renderType << '\n';
        }
        else if (event.key.keysym.sym == SDLK_SPACE) {
            cameraPosition += glm::vec3(0, 0.1, 0)*cameraOrientation;
        }
        else if (event.key.keysym.sym == SDLK_LCTRL) {
            cameraPosition += glm::vec3(0, -0.1, 0)*cameraOrientation;
        }
        else if (event.key.keysym.sym == SDLK_UP) {
            cameraOrientation = cameraOrientation*glm::mat3(
                    1, 0, 0,
                    0, cos(rotationAngle), sin(rotationAngle),
                    0, -sin(rotationAngle), cos(rotationAngle)
            );
        }
        else if (event.key.keysym.sym == SDLK_DOWN) {
            cameraOrientation = cameraOrientation*glm::mat3(
                    1, 0, 0,
                    0, cos(rotationAngle), -sin(rotationAngle),
                    0, sin(rotationAngle), cos(rotationAngle)
            );
        }
        else if (event.key.keysym.sym == SDLK_q) {
            cameraOrientation = cameraOrientation*glm::mat3(
                    cos(rotationAngle), 0, -sin(rotationAngle),
                    0, 1, 0,
                    sin(rotationAngle), 0,  cos(rotationAngle)
            );
        }
        else if (event.key.keysym.sym == SDLK_e) {
            cameraOrientation = cameraOrientation*glm::mat3(
                    cos(rotationAngle), 0, sin(rotationAngle),
                    0, 1, 0,
                    -sin(rotationAngle), 0,  cos(rotationAngle)
            );
        }
    }
    else if (event.key.keysym.sym == SDLK_o) {
        orbit = true;
        std::cout << "Orbit: On" << '\n';
    }
    else if (event.key.keysym.sym == SDLK_l) {
        orbit = false;
        std::cout << "Orbit: Off" << '\n';
    }
    else if (event.key.keysym.sym == SDLK_p) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
    else if (event.type == SDL_MOUSEBUTTONDOWN) {
        std::cout << "Click!" << '\n';
        RayTriangleIntersection rayData;
        RayTriangleIntersection lightData;
        glm::vec3 directionVector = getDirectionVector(glm::vec2(event.button.x, event.button.y), focalLength)*cameraOrientation;
        std::cout << "Coordinates: (" << event.button.x << ", " << event.button.y << ")" << '\n';
        evaluateRay(faces, objects, lights, textures, focalLength, cameraPosition, cameraPosition, directionVector, 1, true);
    }
}

int main(int argc, char *argv[]) {
    OBJSpecs CornellSpecs;
    CornellSpecs.filename = "cornell-box.obj";
    CornellSpecs.clockwise = true;
    CornellSpecs.scaleFactor = 0.17;
    CornellSpecs.offset = glm::vec3(0,0,0);

    OBJSpecs SphereSpecs;
    SphereSpecs.filename = "sphere.obj";
    SphereSpecs.clockwise = true;
    SphereSpecs.scaleFactor = 0.15;
    SphereSpecs.offset = glm::vec3(0,-0.2,0);

    std::string renderType = "raytracing";
    glm::vec3 lightPoint(0,0.4,-0.3); // lights cornell box
//    glm::vec3 lightPoint(0.0,0.5,0.5); // lights top front of sphere
//    glm::vec3 lightPoint(0,0.25,0.2); // lights front of sphere
//    glm::vec3 lightPoint(0.2,0.25,0.4); // lights side of sphere

    glm::vec3 cameraPosition(0, 0, 4);
    int lightNum = 1;
    float lightsize = 0.2; // length of one side of light square
    float focalLength = 2.0;

    glm::mat3 cameraOrientation(
            1,0,0,
            0,1,0,
            0,0,1
            );
    bool orbit = false;

    std::vector<glm::vec3> lights;
    if (std::floor(sqrt(lightNum)) > 1){
        for (float zOffset = -(lightsize/2); zOffset <= lightsize/2; zOffset += lightsize/(std::floor(sqrt(lightNum))-1)) {
            for (float xOffset = -(lightsize/2); xOffset <= lightsize/2; xOffset += lightsize/(std::floor(sqrt(lightNum))-1)) {
                lights.push_back(glm::vec3(lightPoint[0]+xOffset, lightPoint[1], lightPoint[2]+zOffset));
            }
        }
    }
    else {
        lights.push_back(lightPoint);
    }
    std::vector<ModelTriangle> faces;
    std::vector<ModelObject> objects;
    std::map<std::string, TextureMap> textures;

    std::vector<std::vector<float>> depthBuffer = {};
    std::vector<float> column = {};
    for (int i = 0; i < HEIGHT; i++) {
        column.push_back(0.0);
    }
    for (int y = 0; y < WIDTH; y++) {
        depthBuffer.push_back(column);
    }

    importOBJ(CornellSpecs, faces, objects, textures);
//    importOBJ(SphereSpecs, faces, objects, textures);
//    importOBJ(PoolSpecs, faces, objects, textures);
//    importOBJ(LogoSpecs, faces, objects, textures);

    std::map< std::string, std::vector< std::vector< Colour > > > textureGrids;
    for (auto const& texture : textures) {
        textureGrids[texture.first] = {};
        for (int x = 0; x < texture.second.width; x++) {
            textureGrids[texture.first].push_back(std::vector<Colour>());
            for (int y = 0; y < texture.second.height; y++) {
                Colour colour(((texture.second.pixels[((y*texture.second.width)+x)] >> 16) & 0xff), ((texture.second.pixels[((y*texture.second.width)+x)] >> 8) & 0xff), ((texture.second.pixels[((y*texture.second.width)+x)]) & 0xff));
                colour.alpha = ((texture.second.pixels[((y*texture.second.width)+x)] >> 24) & 0xff);
                textureGrids[texture.first][x].push_back(colour);
            }
        }
    }
//    for (int i = 0; i < faces.size(); i++) {
//        std::cout << faces[i] << '\n';
//    }

//    for (int i = 0; i < imagePlane.size(); i++) {
//        std::cout << imagePlane[i] << '\n';
//    }
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
    SDL_Event event;
//    std::vector<glm::vec3> waypoints{glm::vec3(0, 1, 2), glm::vec3(1.5, 1, 0), glm::vec3(0, 2, -2), glm::vec3(-1.5, 3, 0), glm::vec3(0, 1, 2)};
    std::vector<glm::vec3> waypoints{glm::vec3(0, 0, 1), glm::vec3(1, 0, 0), glm::vec3(0, 0, -1), glm::vec3(-1, 0, 0), glm::vec3(0, 0, 1)};
//    std::vector<glm::vec3> waypoints{glm::vec3(0, 0, 4), glm::vec3(0, 0, 1)};
//    std::vector<glm::vec3> waypoints{glm::vec3(0, 0.25, 2), glm::vec3(0.4, 0.35, 1), glm::vec3(-0.4, 0.35, 1), glm::vec3(-0.35, -0.25, 0.5)};
    std::vector<glm::vec3> track;
    std::vector<glm::vec3> trackPart;
    for (int i = 0; i < waypoints.size()-1; i++) {
        trackPart = interpolateThreeElementValues(waypoints[i], waypoints[i+1], 15);
        track.insert(track.end(), trackPart.begin(), trackPart.end());
    }

//    std::cout << "Done!" << '\n';
    while (true) {
        // We MUST poll for events - otherwise the window will freeze !
        if (window.pollForInputEvents(event)) handleEvent(event, window, faces, objects, lights, textureGrids, focalLength, cameraPosition, cameraOrientation, orbit, renderType);
        if (renderType.compare("raytracing") == 0) {
            drawRayTracing(window, faces, objects, lights, textureGrids, focalLength, cameraPosition, cameraOrientation, orbit);
        }
        else if (renderType.compare("rasterise") == 0) {
            drawRasterisedScene(window, faces, textureGrids, depthBuffer, focalLength, cameraPosition, cameraOrientation, orbit);
        }
        else if (renderType.compare("wireframe") == 0) {
            drawWireframe(window, faces, focalLength, cameraPosition, cameraOrientation, orbit);
        };
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }
//    for (int i = 0; i < track.size(); i++) {
////        if (window.pollForInputEvents(event)) handleEvent(event, window, faces, objects, lights, textureGrids, focalLength, cameraPosition, cameraOrientation, orbit, renderType);
//        lookAt(glm::vec3(0,0,0),track[i], cameraOrientation);
////        std::cout << "Camera Position: ("<< track[i][0] << ", " << track[i][1] << ", " << track[i][2] << ") \n";
////        std::cout << "(" << cameraOrientation[0][0] << ", " << cameraOrientation[1][0] << ", " << cameraOrientation[2][0] <<  '\n';
////        std::cout << "," << cameraOrientation[0][1] << ", " << cameraOrientation[1][1] << ", " << cameraOrientation[2][1] <<  '\n';
////        std::cout << "," << cameraOrientation[0][2] << ", " << cameraOrientation[1][2] << ", " << cameraOrientation[2][2] << ")" <<  '\n' <<  '\n';
////        if (i < (50)){
////            drawWireframe(window, faces, focalLength, track[i], cameraOrientation, false);
////        }
////        else if (i < 100) {
////            drawRasterisedScene(window, faces, textureGrids, depthBuffer, focalLength, track[i], cameraOrientation, false);
////        }
////        else {
//            drawRayTracing(window, faces, objects, lights, textureGrids, focalLength, track[i], cameraOrientation, false);
////        }
//        window.renderFrame();
//        window.saveBMP("output/"+std::to_string(i)+".bmp");
//    }
    return 0;
}