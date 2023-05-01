#include "hw.hpp"

#include <iostream>
#include <vector>

namespace COL781 {
	namespace OpenGL {

		GLenum glCheckError_(const char *file, int line) {
			GLenum errorCode;
			while ((errorCode = glGetError()) != GL_NO_ERROR) {
				std::string error;
				switch (errorCode) {
				case GL_INVALID_ENUM:
					error = "INVALID_ENUM";
					break;
				case GL_INVALID_VALUE:
					error = "INVALID_VALUE";
					break;
				case GL_INVALID_OPERATION:
					error = "INVALID_OPERATION";
					break;
				case GL_STACK_OVERFLOW:
					error = "STACK_OVERFLOW";
					break;
				case GL_STACK_UNDERFLOW:
					error = "STACK_UNDERFLOW";
					break;
				case GL_OUT_OF_MEMORY:
					error = "OUT_OF_MEMORY";
					break;
				case GL_INVALID_FRAMEBUFFER_OPERATION:
					error = "INVALID_FRAMEBUFFER_OPERATION";
					break;
				}
				std::cout << error << " | " << file << " (" << line << ")" << std::endl;
			}
			return errorCode;
		}
#define glCheckError() glCheckError_(__FILE__, __LINE__) 

		bool Rasterizer::initialize(const std::string &title, int width, int height, int spp) {
			if (SDL_Init(SDL_INIT_VIDEO) < 0) {
				std::cout << "Could not initialize SDL: " << SDL_GetError() << std::endl;
				return false;
			}
			SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
			SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
			SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
			SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
			SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, spp);
			window = SDL_CreateWindow(title.c_str(), SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height, SDL_WINDOW_OPENGL);
			if (!window) {
				std::cerr << "Could not create window: " << SDL_GetError() << std::endl;
				return false;
			}
			if (!SDL_GL_CreateContext(window)) {
				std::cerr << "Could not create OpenGL context: " << SDL_GetError() << std::endl;
				return false;
			}
			GLenum glewStatus = glewInit();
			if (glewStatus != GLEW_OK) {
				std::cerr << "Could not initialize GLEW: " << glewGetErrorString(glewStatus) << std::endl;
				return false;
			}
			quit = false;
			glCheckError();
			return true;
		}

		bool Rasterizer::shouldQuit() {
			glCheckError();
			return quit;
		}

		ShaderProgram Rasterizer::createShaderProgram(const VertexShader &vs, const FragmentShader &fs) {
			ShaderProgram program = glCreateProgram();
			glAttachShader(program, vs);
			glAttachShader(program, fs);
			glLinkProgram(program);
			GLint linkStatus;
			glGetProgramiv(program, GL_LINK_STATUS, &linkStatus);
			if (linkStatus != GL_TRUE) {
				std::cout << "Error linking shaders:" << std::endl;
				GLint maxLength = 0;
				glGetProgramiv(program, GL_INFO_LOG_LENGTH, &maxLength);
				std::vector<GLchar> infoLog(maxLength);
				glGetShaderInfoLog(fs, maxLength, &maxLength, &infoLog[0]);
				std::cout << &infoLog[0] << std::endl;
				glDeleteProgram(program);
				glDeleteShader(vs);
				glDeleteShader(fs);
				return 0;
			}
			glCheckError();
			return program;
		}

		void Rasterizer::useShaderProgram(const ShaderProgram &program) {
			glUseProgram(program);
			glCheckError();
		}

		template <> void Rasterizer::setUniform(ShaderProgram &program, const std::string &name, float value) {
			GLint location = glGetUniformLocation(program, name.c_str());
			glUniform1f(location, value);
			glCheckError();
		}
		
		template <> void Rasterizer::setUniform(ShaderProgram &program, const std::string &name, int value) {
			GLint location = glGetUniformLocation(program, name.c_str());
			glUniform1i(location, value);
			glCheckError();
		}

		template <> void Rasterizer::setUniform(ShaderProgram &program, const std::string &name, glm::vec2 value) {
			GLint location = glGetUniformLocation(program, name.c_str());
			glUniform2fv(location, 1, &value[0]);
			glCheckError();
		}
		
		template <> void Rasterizer::setUniform(ShaderProgram &program, const std::string &name, glm::vec3 value) {
			GLint location = glGetUniformLocation(program, name.c_str());
			glUniform3fv(location, 1, &value[0]);
			glCheckError();
		}
		
		template <> void Rasterizer::setUniform(ShaderProgram &program, const std::string &name, glm::vec4 value) {
			GLint location = glGetUniformLocation(program, name.c_str());
			glUniform4fv(location, 1, &value[0]);
			glCheckError();
		}

		template <> void Rasterizer::setUniform(ShaderProgram &program, const std::string &name, glm::mat2 value) {
			GLint location = glGetUniformLocation(program, name.c_str());
			glUniformMatrix2fv(location, 1, GL_FALSE, &value[0][0]);
			glCheckError();
		}

		template <> void Rasterizer::setUniform(ShaderProgram &program, const std::string &name, glm::mat3 value) {
			GLint location = glGetUniformLocation(program, name.c_str());
			glUniformMatrix3fv(location, 1, GL_FALSE, &value[0][0]);
			glCheckError();
		}

		template <> void Rasterizer::setUniform(ShaderProgram &program, const std::string &name, glm::mat4 value) {
			GLint location = glGetUniformLocation(program, name.c_str());
			glUniformMatrix4fv(location, 1, GL_FALSE, &value[0][0]);
			glCheckError();
		}

		void Rasterizer::deleteShaderProgram(ShaderProgram &program) {
			glDeleteProgram(program);
			glCheckError();
		}

		Object Rasterizer::createObject() {
			Object object;
			glGenVertexArrays(1, &object.vao);
			glCheckError();
			return object;
		}

		// template <typename T> int dim();
		// template <> int dim<float>() {return 1;}
		// template <> int dim<glm::vec2>() {return 2;}
		// template <> int dim<glm::vec3>() {return 3;}
		// template <> int dim<glm::vec4>() {return 4;}

		template <typename T> AttribBuf Rasterizer::createVertexAttribs(Object &object, int attribIndex, int n, const T *data) {
			GLuint vbo;
			glGenBuffers(1, &vbo);
			glBindVertexArray(object.vao);
			glBindBuffer(GL_ARRAY_BUFFER, vbo);
			glBufferData(GL_ARRAY_BUFFER, n*sizeof(T), (float*)data, GL_STATIC_DRAW);
			glVertexAttribPointer(attribIndex, sizeof(T)/sizeof(float), GL_FLOAT, GL_FALSE, sizeof(T), NULL);
			glEnableVertexAttribArray(attribIndex);
			glCheckError();
			return vbo;
		}
		template AttribBuf Rasterizer::createVertexAttribs(Object &object, int attribIndex, int n, const float *data);
		template AttribBuf Rasterizer::createVertexAttribs(Object &object, int attribIndex, int n, const glm::vec2 *data);
		template AttribBuf Rasterizer::createVertexAttribs(Object &object, int attribIndex, int n, const glm::vec3 *data);
		template AttribBuf Rasterizer::createVertexAttribs(Object &object, int attribIndex, int n, const glm::vec4 *data);

		template <typename T> void Rasterizer::updateVertexAttribs(AttribBuf &vbo, int n, const T* data) {
			glBindBuffer(GL_ARRAY_BUFFER, vbo);
			glBufferData(GL_ARRAY_BUFFER, n*sizeof(T), (float*)data, GL_STATIC_DRAW);
			glCheckError();
		}
		template void Rasterizer::updateVertexAttribs(AttribBuf &vbo, int n, const float* data);
		template void Rasterizer::updateVertexAttribs(AttribBuf &vbo, int n, const glm::vec2* data);
		template void Rasterizer::updateVertexAttribs(AttribBuf &vbo, int n, const glm::vec3* data);
		template void Rasterizer::updateVertexAttribs(AttribBuf &vbo, int n, const glm::vec4* data);

		IndexBuf Rasterizer::createTriangleIndices(Object &object, int n, const glm::ivec3* indices) {
			GLuint ebo;
			glGenBuffers(1, &ebo);
			glBindVertexArray(object.vao);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*n*sizeof(int), (float*)indices, GL_STATIC_DRAW);
			object.nTris = n;
			glCheckError();
			return ebo;
		}
		
		void Rasterizer::enableDepthTest() {
			glEnable(GL_DEPTH_TEST);
		    glDepthFunc(GL_LESS);   
			glCheckError();
		}

		void Rasterizer::clear(glm::vec4 color) {
			glClearColor(color[0], color[1], color[2], color[3]);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glCheckError();
		}

		// template <> Buffer<glm::vec4> Rasterizer::bufferVertexData(int n, glm::vec4* data){
		// 	GLuint VBO;
		// 	glGenBuffers(1, &VBO);
		// 	glBindBuffer(GL_ARRAY_BUFFER, VBO);
		// 	glBufferData(GL_ARRAY_BUFFER, n*sizeof(glm::vec4), data, GL_STATIC_DRAW);
		// 	return VBO;
		// }

		// Buffer<glm::ivec3> Rasterizer::bufferElements(int n, glm::ivec3* data) {
		// 	GLuint IBO;
		// 	glGenBuffers(1, &IBO);
		// 	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
		// 	glBufferData(GL_ELEMENT_ARRAY_BUFFER, n * sizeof(glm::ivec3), data, GL_STATIC_DRAW);
		// 	return IBO;
		// }

		// template <> void Rasterizer::setVertexAttribs<glm::vec4>(const ShaderProgram &program, const std::string &name, const Buffer<glm::vec4> buffer) {

		// 	GLint loc = glGetAttribLocation(program, name.c_str());

		// 	glEnableVertexAttribArray(loc);

		// 	glBindBuffer(GL_ARRAY_BUFFER, buffer);
		// 	glVertexAttribPointer(
		// 						  loc,
		// 						  4,
		// 						  GL_FLOAT,
		// 						  GL_FALSE,
		// 						  0, 
		// 						  (void *)0
		// 						  );

		// }

		void Rasterizer::drawObject(const Object &object) {
			glBindVertexArray(object.vao);
			glDrawElements(GL_TRIANGLES, 3*object.nTris, GL_UNSIGNED_INT, 0);
			glCheckError();
		}

		void Rasterizer::setupFilledFaces() {
	        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}

		void Rasterizer::setupWireFrame() {
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glEnable(GL_POLYGON_OFFSET_LINE);
			glPolygonOffset(-1.f, -1.f);
		}

		void Rasterizer::show() {
			SDL_GL_SwapWindow(window);
			SDL_Event e;
			while (SDL_PollEvent(&e) != 0) {
				if(e.type == SDL_QUIT) {
					quit = true;
				}
			}
			glCheckError();
		}

		// glm::vec3 Rasterizer::getCameraUpdate(float cameraSpeed) {
		// 	SDL_Event e;
		// 	while (SDL_PollEvent(&e) != 0) {
		// 		if(e.type == SDL_KEYDOWN) {
		// 			std::cout<<e.key.keysym.sym<<" key pressed"<<std::endl;
		// 		}
		// 	}
		// 	return glm::vec3(0.0f);
		// }
		// glm::vec3 Rasterizer::getCameraUpdate() {
		// 	SDL_GL_SwapWindow(window);
		// 	SDL_Event e;
		// 	int initialPosx, initialPosy, finalPosx, finalPosy;
		// 	bool camDragEvent = false;

		// 	while (SDL_PollEvent(&e) != 0) {

		// 		if(e.type==SDL_MOUSEBUTTONDOWN) {
		// 			SDL_GetGlobalMouseState(&initialPosx, &initialPosy);
		// 			std::cout<<"Mouse click "<<" xM "<<initialPosx<<" yM "<<initialPosy<<std::endl;
		// 			while( e.type!=SDL_MOUSEBUTTONUP);
					
		// 				SDL_GetGlobalMouseState(&finalPosx, &finalPosy);
		// 				std::cout<<"Mouse release "<<" xM "<<finalPosx<<" yM "<<finalPosy<<std::endl;
		// 				camDragEvent = true;

		// 		}
		// 	}
		// 	glCheckError();
		// 	if(camDragEvent) {
		// 		return 0.001f *glm::vec3((float)(finalPosx - initialPosx), (float)(finalPosy - initialPosy), 0.0f);
		// 	}
		// 	else {
		// 		return glm::vec3(0.0f, 0.0f, 0.0f);
		// 	}
		// }

		// TODO WIP
		/*glm::vec3 Rasterizer::getCameraUpdate(float x, float y) {
			float lastX = 640 / 2, lastY = 720 / 2;
			float pitch = 0.0f, yaw = -90.0f;

			float offsetX = x - lastX;
			float offsetY = lastY - y;
			lastX = x;
			lastY = y;

			float sensitivity = 0.3f;

			yaw += offsetX * sensitivity;
			pitch += offsetY * sensitivity;

			if (pitch > 89.0f)
				pitch = 89.0f;
			if (pitch < -89.0f)
				pitch = -89.0f;

			glm::vec3 front;
			front.x = std::cos(glm::radians(yaw)) * std::cos(glm::radians(pitch));
			front.y = std::sin(glm::radians(pitch));
			front.z = std::sin(glm::radians(yaw)) * std::cos(glm::radians(pitch));

			//std::cout << "Mouse: " << offsetX << ", " << offsetY << std::endl;
			//std::cout << "YAW AND PITCH: " << yaw << ", " << pitch << std::endl;

			return glm::normalize(front);
		}*/

		GLuint createShader(GLenum type, const char *source) {
			GLuint shader = glCreateShader(type);
			glShaderSource(shader, 1, &source, NULL);
			glCompileShader(shader);
			GLint compileStatus;
			glGetShaderiv(shader, GL_COMPILE_STATUS, &compileStatus);
			if (compileStatus != GL_TRUE) {
				std::cout << "Error compiling shader:" << std::endl;
				GLint maxLength = 0;
				glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &maxLength);
				std::vector<GLchar> infoLog(maxLength);
				glGetShaderInfoLog(shader, maxLength, &maxLength, &infoLog[0]);
				std::cout << &infoLog[0] << std::endl;
				glDeleteShader(shader);
				return 0;
			}
			glCheckError();
			return shader;
		}

		VertexShader Rasterizer::vsPhongShading() {
			const char *source =
				"#version 330 core\n"
				"layout(location = 0) in vec3 vertex;\n"
				"layout(location = 1) in vec3 normal;\n"
				"uniform mat4 model;\n"
				"uniform mat4 view;\n"
				"uniform mat4 projection;\n"
				"out vec3 FragPos;\n"
				"out vec3 Normal;\n"
				"void main() {\n"
				"FragPos = mat3(model) * vertex;\n"			
				"Normal = transpose(inverse(mat3(model))) * normal;\n"
				"gl_Position = projection * view * model * vec4(vertex,1.0);\n"
				"}\n";

			return createShader(GL_VERTEX_SHADER, source);
		}

		FragmentShader Rasterizer::fsPhongShading() {
			const char *source =
				"#version 330 core\n"  
				"in vec3 FragPos;\n"
				"in vec3 Normal;\n"
				"out vec4 fColor;\n"
				"uniform vec3 lightPos;\n"
				"uniform vec3 viewPos;\n"
				"uniform vec3 lightColor;\n"
				"uniform vec3 objectColor;\n"
				"void main() {\n"
				"// ambient\n"
				"float Ka = 0.4;\n"
				"vec3 ambient = vec3(Ka);\n"
				"// diffuse\n"
				"float Kd = 0.5;\n"
				"vec3 norm = normalize(Normal);\n"
				"vec3 lightDir = normalize(lightPos - FragPos);\n"
				"float diff = max(dot(norm, lightDir), 0.0); \n"
				"vec3 diffuse = Kd * diff * lightColor;\n"
				"// specular \n"
				"float Ks = 0.1;\n"
				"float p = 64;\n"
				"vec3 viewDir = normalize(viewPos - FragPos);\n"
				"vec3 halfDir = normalize(viewDir + lightDir); \n"
				// "vec3 reflectDir = reflect(-lightDir, norm); \n"
				"float spec = pow(max(dot(halfDir, norm), 0.0), p);\n"
				"vec3 specular = Ks * spec * lightColor;\n"
				"vec3 result = (ambient + diffuse) * objectColor + specular; \n"
				"fColor = vec4(result, 1.0);\n"
				"}\n";
			return createShader(GL_FRAGMENT_SHADER, source);
		}
			
	}
}
