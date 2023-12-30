//
// Created by Pawulon on 30/12/2023.
//

#include "CreateOpenGlWindow.hpp"
//#ifndef STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "glad/stb_image_write.h"
//#endif

void saveScreenshot(const char* filename, GLFWwindow* window) {

    int width, height;

    glfwGetWindowSize(window, &width, &height);
    unsigned char* pixels = new unsigned char[3 * width * height];

    // Read the pixel data from the framebuffer
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

    // Flip the image vertically
    for (int i = 0; i < height / 2; ++i) {
        for (int j = 0; j < width * 3; ++j) {
            std::swap(pixels[i * width * 3 + j], pixels[(height - i - 1) * width * 3 + j]);
        }
    }

    // Save the image using stb_image_write
    stbi_write_jpg(filename, width, height, 3, pixels, 100);

    delete[] pixels;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

int CreateOpenGlWindow(std::vector<Vertex2D> &vertices, std::vector<fem::TriangleElement> &Elements, std::vector<double> &c, double L){
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    GLFWwindow* window = glfwCreateWindow(800, 800, "LearnOpenGL", nullptr, nullptr);
    if (window == nullptr)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }


    auto shaderProgram = createShaderProgram();


    // Vertex Data
    std::vector<GLfloat> verticesBuffer;

    std::vector<GLuint> indecesBuffer;



    for(int i=0; i<vertices.size(); i++){
        verticesBuffer.insert(verticesBuffer.end(), {static_cast<GLfloat>(vertices[i].x),
                                                     static_cast<GLfloat>(vertices[i].y),
                                                     static_cast<GLfloat>(c[i])});

    }

    for(auto&e:Elements){
        for(auto &i:e.globalVectorIdx){
            indecesBuffer.insert(indecesBuffer.end(), {static_cast<GLuint>(i)});
        }
    }

//    vertices.insert(vertices.end(), {0.5f, -0.5f, 0.0f, 1.0f,  0.5f, 0.5f, 0.0f,  0.5f, 0.5f});

    // Vertex Buffer Object (VBO)
    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    // Bind VAO
    glBindVertexArray(VAO);

    // Bind VBO and copy vertices data into it
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, verticesBuffer.size() * sizeof(GLfloat) , &verticesBuffer[0], GL_STATIC_DRAW);

    //
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indecesBuffer.size() * sizeof(GLuint) , &indecesBuffer[0], GL_STATIC_DRAW);

    // Set vertex attribute pointers
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), NULL);
    glEnableVertexAttribArray(0);

    // Unbind VAO and VBO
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glViewport(0, 0, 800, 800);
    glClearColor(0.1, 0.2, 0.2, 0.0);


    auto ortMtrx = makeOrto(-L/2, L/2,-L/2, L/2, -2., 12.);


    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    while(!glfwWindowShouldClose(window))
    {
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, true);

        if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS)
        {
            saveScreenshot("screenshot.jpg", window);

        }


        // Rendering commands
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // Use the shader program
        glUseProgram(shaderProgram);

        int projectionLocation = glGetUniformLocation(shaderProgram, "projection");

        glUniformMatrix4fv(projectionLocation, 1, GL_FALSE, ortMtrx.get());

        // Bind the VAO
        glBindVertexArray(VAO);

        // Draw the triangle
        glDrawElements(GL_TRIANGLES, indecesBuffer.size(), GL_UNSIGNED_INT, 0);

        // Unbind the VAO
        glBindVertexArray(0);

        // Swap buffers and poll events
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Clean up resources
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shaderProgram);

    glfwTerminate();
    return 0;
}
