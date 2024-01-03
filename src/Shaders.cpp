//
// Created by Pawulon on 30/12/2023.
//

#include "../Visualization/Shaders.hpp"

unsigned int createShaderProgram(){

    const char* vertexShaderSource = R"(
        #version 330 core
        in  vec3 position;

        uniform mat4 projection;

        out float zPosition;

        void main() {
            gl_Position =  projection * vec4(position, 1.0);
            zPosition = position.z;
        }
    )";

    // Fragment Shader Source Code
    const char* fragmentShaderSource = R"(
        #version 330 core

        in float zPosition;
        out vec4 fragColor;

        void main() {
            float t = zPosition;
            vec3 color1 = vec3(0.0, 1.0, 1.0);
            vec3 color2 = vec3(0.0, 0.0, 0.0);
            vec3 color3 = vec3(1., 1.0, 1.);

//            vec3 color = mix(color2, color3 , log(t + 1)/log(2));
            vec3 color = mix(color2, color3 , t);
//            if(t < 0.5){
//                 color = mix(color1, color2 , t/.5);
//            }else{
//                 color = mix(color2, color3 , (t -.5)*2);
//            }

            fragColor = vec4(color, 1.0);
        }
    )";

    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
    glCompileShader(vertexShader);

    int success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertexShader, 512, nullptr, infoLog);
        std::cerr << "Vertex shader compilation failed:\n" << infoLog << std::endl;
    }

    // Fragment Shader
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
    glCompileShader(fragmentShader);

    // Check for compile errors
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragmentShader, 512, nullptr, infoLog);
        std::cerr << "Fragment shader compilation failed:\n" << infoLog << std::endl;
    }

    // Shader Program
    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    // Check for linking errors
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, nullptr, infoLog);
        std::cerr << "Shader program linking failed:\n" << infoLog << std::endl;
    }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    return shaderProgram;
}

std::shared_ptr<GLfloat> makeOrto(double l, double r, double b, double t, double f, double n){
    auto* otrh_mtx = new GLfloat[16]{
            static_cast<GLfloat>(2./(r-l)),            static_cast<GLfloat>(0),                static_cast<GLfloat>(0),                static_cast<GLfloat>(0),
            static_cast<GLfloat>(0),                   static_cast<GLfloat>(2./(t-b)),         static_cast<GLfloat>(0),                static_cast<GLfloat>(0),
            static_cast<GLfloat>(0),                   static_cast<GLfloat>(0),                static_cast<GLfloat>(-2./(f-n)),        static_cast<GLfloat>(0),
            static_cast<GLfloat>(-(r+l)/(r-l)),        static_cast<GLfloat>(-(t+b)/(t-b)),     static_cast<GLfloat>(-(f+n)/(f-n)),     static_cast<GLfloat>(1)
    };

    return {otrh_mtx, [](const GLfloat* arr) {
        delete[] arr; // Custom deleter for the array
    }};
}