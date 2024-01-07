//
// Created by Pawulon on 30/12/2023.
//

#include "../Visualization/Shaders.hpp"

namespace mes::plot{
unsigned int createShaderProgram(bool log) {

    const char *vertexShaderSource = R"(
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

    const char *fragmentShaderSourceLogColor = R"(
        #version 330 core

        float colormap_f1(float x) {
            return -510.0 * x + 255.0;
        }

        float colormap_f2(float x) {
            return (-1891.7 * x + 217.46) * x + 255.0;
        }

        float colormap_f3(float x) {
            return 9.26643676359015e1 * sin((x - 4.83450094847127e-1) * 9.93) + 1.35940451627965e2;
        }

        float colormap_f4(float x) {
            return -510.0 * x + 510.0;
        }

        float colormap_f5(float x) {
            float xx = x - 197169.0 / 251000.0;
            return (2510.0 * xx - 538.31) * xx;
        }

        float colormap_red(float x) {
            if (x < 0.0) {
                return 1.0;
            } else if (x < 10873.0 / 94585.0) {
                float xx = colormap_f2(x);
                if (xx > 255.0) {
                    return (510.0 - xx) / 255.0;
                } else {
                    return xx / 255.0;
                }
            } else if (x < 0.5) {
                return 1.0;
            } else if (x < 146169.0 / 251000.0) {
                return colormap_f4(x) / 255.0;
            } else if (x < 197169.0 / 251000.0) {
                return colormap_f5(x) / 255.0;
            } else {
                return 0.0;
            }
        }

        float colormap_green(float x) {
            if (x < 10873.0 / 94585.0) {
                return 1.0;
            } else if (x < 36373.0 / 94585.0) {
                return colormap_f2(x) / 255.0;
            } else if (x < 0.5) {
                return colormap_f1(x) / 255.0;
            } else if (x < 197169.0 / 251000.0) {
                return 0.0;
            } else if (x <= 1.0) {
                return abs(colormap_f5(x)) / 255.0;
            } else {
                return 0.0;
            }
        }

        float colormap_blue(float x) {
            if (x < 0.0) {
                return 0.0;
            } else if (x < 36373.0 / 94585.0) {
                return colormap_f1(x) / 255.0;
            } else if (x < 146169.0 / 251000.0) {
                return colormap_f3(x) / 255.0;
            } else if (x <= 1.0) {
                return colormap_f4(x) / 255.0;
            } else {
                return 0.0;
            }
        }

        vec4 colormap(float x) {
            return vec4(colormap_red(x), colormap_green(x), colormap_blue(x), 1.0);
        }

        in float zPosition;
        out vec4 fragColor;

        void main() {

            float value = zPosition;
            value = 1. - log(zPosition * 9. + 1.)/log(10) ;
            //value = log(value*99. + 1.)/log(100.);
            vec3 color1 = vec3(0.0, 0.0, 0.0);
            vec3 color2 = vec3(1., 1., 1.);


            vec3 color = mix(color1, color2 , value);

            fragColor = colormap(value);
        }
    )";

    const char *fragmentShaderSource = R"(
        #version 330 core

        in float zPosition;
        out vec4 fragColor;

        void main() {
            float value = zPosition;
            vec3 color1 = vec3(0.,0.,0.);
            vec3 color2 = vec3(1, 0., 1 );
            vec3 color3 = vec3(0, 1, 1);

            vec3 color = mix(color1, color2 , log(((value -.5)*2.)*9.+ 1.)/log(10.));
            if(value <.5){
                 color = mix(color1, color3 , log(((.5 - value)*2.)*9.+ 1.)/log(10.));
            }

            gl_FragColor = vec4(color, 1.0);
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
    if(log){
        glShaderSource(fragmentShader, 1, &fragmentShaderSourceLogColor, nullptr);
    }else{
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
    }

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

    unsigned int createShaderProgramGeometry(bool log) {


        const char* vertexShaderSource = R"(
            #version 330 core
            in vec3 position;

        //    uniform mat4 projection;

        //    out DATA
        //    {
        //        mat4 ortProjection;
        //    } data_out;

            void main() {
                gl_Position =  vec4(position, 1.0);
        //        data_out.ortProjection = projection;
            }
        )";

        const char* geometryShaderSource = R"(
            #version 330 core

            layout (triangles_adjacency) in;
            layout (triangle_strip, max_vertices = 8) out;

            uniform mat4 projection;

            out float zPosition;

            in DATA
            {
                mat4 ortProjection;
            } data_in[];

            float lin_phi0(float zeta, float eta) {
                return -.5 * (zeta + eta);
            }

            float lin_phi1(float zeta, float eta) {
                return .5 * (1 + zeta);
            }

            float lin_phi2(float zeta, float eta) {
                return .5 * (1 + eta);
            }

            float quad_phi0(float zeta, float eta) {
                float L = lin_phi0(zeta, eta);
                return  L * (2. * L - 1);
            }


            float quad_phi1(float zeta, float eta) {
                float L = lin_phi1(zeta, eta);
                return  L * (2.f * L - 1.f);
            }

            float quad_phi2(float zeta, float eta) {
                float L = lin_phi2(zeta, eta);
                return  L * (2. * L - 1);
            }

            float quad_phi3(float zeta, float eta) {
                return  4. * lin_phi0(zeta, eta)*lin_phi1(zeta, eta);

            }

            float quad_phi4(float zeta, float eta) {
                return  4. * lin_phi1(zeta, eta)*lin_phi2(zeta, eta);

            }
            float quad_phi5(float zeta, float eta) {
                return  4. * lin_phi0(zeta, eta)*lin_phi2(zeta, eta);
            }

            float mapx(float x, float y){
                return        gl_in[0].gl_Position.x*lin_phi0(x, y)
                            + gl_in[1].gl_Position.x*lin_phi1(x, y)
                            + gl_in[2].gl_Position.x*lin_phi2(x, y);
            }

            float mapy(float x, float y){
                return        gl_in[0].gl_Position.y*lin_phi0(x, y)
                            + gl_in[1].gl_Position.y*lin_phi1(x, y)
                            + gl_in[2].gl_Position.y*lin_phi2(x, y);
            }

            float mapu(float x, float y){
                return        gl_in[0].gl_Position.z*quad_phi0(x, y)
                            + gl_in[1].gl_Position.z*quad_phi1(x, y)
                            + gl_in[2].gl_Position.z*quad_phi2(x, y)
                            + gl_in[3].gl_Position.z*quad_phi3(x, y)
                            + gl_in[4].gl_Position.z*quad_phi4(x, y)
                            + gl_in[5].gl_Position.z*quad_phi5(x, y);
            }

            vec2 calcU(float x, float y){
                return vec2(mapx(x,y), mapy(x,y));
            }

            // Default main function
            void main()
            {
                vec2 tempPos = vec2(-1.f, 1.f);


                int n = 2;
                float step = 2.f/n;
                for(int i=0; i<n; i++){
                    tempPos = calcU(-1.f + i*step, 1.f - i*step);
                    gl_Position = projection * vec4(tempPos, 1.f, 1.f);
                    zPosition = mapu(-1.f + i*step, 1.f - i*step);
                    EmitVertex();
                    for(int j = i+1; j <= n; j++){
                        tempPos = calcU(-1.f + i*step, 1.f - j*step);
                        gl_Position = projection * vec4(tempPos, 1.f, 1.f);
                        zPosition = mapu(-1.f + i*step, 1.f - j*step);
                        EmitVertex();

                        tempPos = calcU(-1.f + (i + 1.f) * step, 1.f - j*step);
                        gl_Position = projection * vec4(tempPos, 1.f, 1.f);
                        zPosition = mapu(-1.f + (i + 1.f) * step, 1.f - j*step);
                        EmitVertex();

                    }
                    EndPrimitive();
                }
            }
        )";

        const char *fragmentShaderSource = R"(
        #version 330 core

        float colormap_f1(float x) {
            return -510.0 * x + 255.0;
        }

        float colormap_f2(float x) {
            return (-1891.7 * x + 217.46) * x + 255.0;
        }

        float colormap_f3(float x) {
            return 9.26643676359015e1 * sin((x - 4.83450094847127e-1) * 9.93) + 1.35940451627965e2;
        }

        float colormap_f4(float x) {
            return -510.0 * x + 510.0;
        }

        float colormap_f5(float x) {
            float xx = x - 197169.0 / 251000.0;
            return (2510.0 * xx - 538.31) * xx;
        }

        float colormap_red(float x) {
            if (x < 0.0) {
                return 1.0;
            } else if (x < 10873.0 / 94585.0) {
                float xx = colormap_f2(x);
                if (xx > 255.0) {
                    return (510.0 - xx) / 255.0;
                } else {
                    return xx / 255.0;
                }
            } else if (x < 0.5) {
                return 1.0;
            } else if (x < 146169.0 / 251000.0) {
                return colormap_f4(x) / 255.0;
            } else if (x < 197169.0 / 251000.0) {
                return colormap_f5(x) / 255.0;
            } else {
                return 0.0;
            }
        }

        float colormap_green(float x) {
            if (x < 10873.0 / 94585.0) {
                return 1.0;
            } else if (x < 36373.0 / 94585.0) {
                return colormap_f2(x) / 255.0;
            } else if (x < 0.5) {
                return colormap_f1(x) / 255.0;
            } else if (x < 197169.0 / 251000.0) {
                return 0.0;
            } else if (x <= 1.0) {
                return abs(colormap_f5(x)) / 255.0;
            } else {
                return 0.0;
            }
        }

        float colormap_blue(float x) {
            if (x < 0.0) {
                return 0.0;
            } else if (x < 36373.0 / 94585.0) {
                return colormap_f1(x) / 255.0;
            } else if (x < 146169.0 / 251000.0) {
                return colormap_f3(x) / 255.0;
            } else if (x <= 1.0) {
                return colormap_f4(x) / 255.0;
            } else {
                return 0.0;
            }
        }

        vec4 colormap(float x) {
            return vec4(colormap_red(x), colormap_green(x), colormap_blue(x), 1.0);
        }

        in float zPosition;
        out vec4 fragColor;

        void main() {

            float value = zPosition;
            value = 1. - log(zPosition * 9. + 1.)/log(10) ;
            //value = log(value*99. + 1.)/log(100.);
            vec3 color1 = vec3(0.0, 0.0, 0.0);
            vec3 color2 = vec3(1., 1., 1.);


            vec3 color = mix(color1, color2 , value);

            fragColor = colormap(value);
        }
    )";

        const char *fragmentShaderSource2 = R"(
        #version 330 core

        in float zPosition;
        out vec4 fragColor;


        void main() {
            float value = zPosition;
            vec3 color1 = vec3(0.,0.,0.);
            vec3 color2 = vec3(1, 0., 1 );
            vec3 color3 = vec3(0, 1, 1);

            vec3 color = mix(color1, color2 , log(((value -.5)*2.)*9.+ 1.)/log(10.));
            if(value <.5){
                 color = mix(color1, color3 , log(((.5 - value)*2.)*9.+ 1.)/log(10.));
            }

            gl_FragColor = vec4(color, 1.0);
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
//        glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
        if(log){
            glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
        }else{
            glShaderSource(fragmentShader, 1, &fragmentShaderSource2, nullptr);
        }
        glCompileShader(fragmentShader);

        // Check for compile errors
        glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
        if (!success) {
            glGetShaderInfoLog(fragmentShader, 512, nullptr, infoLog);
            std::cerr << "Fragment shader compilation failed:\n" << infoLog << std::endl;
        }

        unsigned int geometryShader = glCreateShader(GL_GEOMETRY_SHADER);
        // Attach Geometry Shader source to the Fragment Shader Object
        glShaderSource(geometryShader, 1, &geometryShaderSource, nullptr);
        // Compile the Geometry Shader into machine code
        glCompileShader(geometryShader);

        glGetShaderiv(geometryShader, GL_COMPILE_STATUS, &success);
        if (!success) {
            glGetShaderInfoLog(geometryShader, 512, nullptr, infoLog);
            std::cerr << "Geometry shader compilation failed:\n" << infoLog << std::endl;
        }

        // Shader Program
        unsigned int shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        glAttachShader(shaderProgram, geometryShader);
        glLinkProgram(shaderProgram);

        // Check for linking errors
        glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
        if (!success) {
            glGetProgramInfoLog(shaderProgram, 512, nullptr, infoLog);
            std::cerr << "Shader program linking failed:\n" << infoLog << std::endl;
        }

        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
        glDeleteShader(geometryShader);

        return shaderProgram;
    }


std::shared_ptr<GLfloat> makeOrto(double l, double r, double b, double t, double f, double n) {
    auto *otrh_mtx = new GLfloat[16]{
            static_cast<GLfloat>(2. / (r - l)), static_cast<GLfloat>(0), static_cast<GLfloat>(0),
            static_cast<GLfloat>(0),
            static_cast<GLfloat>(0), static_cast<GLfloat>(2. / (t - b)), static_cast<GLfloat>(0),
            static_cast<GLfloat>(0),
            static_cast<GLfloat>(0), static_cast<GLfloat>(0), static_cast<GLfloat>(-2. / (f - n)),
            static_cast<GLfloat>(0),
            static_cast<GLfloat>(-(r + l) / (r - l)), static_cast<GLfloat>(-(t + b) / (t - b)),
            static_cast<GLfloat>(-(f + n) / (f - n)), static_cast<GLfloat>(1)
    };

    return {otrh_mtx,
            [](const GLfloat *arr) {
                delete[] arr; // Custom deleter for the array
            }
    };
}

}