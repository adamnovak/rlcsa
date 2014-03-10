package fi.helsinki.cs.rlcsa;

import java.io.*;
import java.nio.file.Files;

/**
 * A class to load the RLCSA native library whenever the RLCSA jar is loaded.
 * Only works if it's packaged into a jar with "rlcsa.so" in it at the root
 * level.
 *
 * Based on:
 * Simple library class for working with JNI (Java Native Interface)
 * @see http://frommyplayground.com/how-to-load-native-jni-library-from-jar
 * @author Adam Heirnich <adam@adamh.cz>, http://www.adamh.cz
 * Licensed under the MIT license.
 */
class RLCSANativeLoader {
    static {
        try {
            // Read the library data out of our jar
        	InputStream libraryData = 
        	    RLCSANativeLoader.class.getResourceAsStream("rlcsa.so");
            if (libraryData == null) {
                // Complain we didn't find the library in our jar. See
                // <http://frommyplayground.com/how-to-load-native-jni-library-
                // from-jar/>
                throw new FileNotFoundException(
                    "RLCSA native library not in jar!");
            }
            
            // Make a temporary file to hold the library
            File library = File.createTempFile("rlcsa", "so");
            // Mark it to be deleted when we exit.
            library.deleteOnExit();
            
            // Copy all the data over with the new NIO file copying method.
            Files.copy(libraryData, library.toPath());
            
            // Load the library
        	System.load(library.getAbsolutePath());
	    } catch(Exception e) {
	        // Catch any ordinary "checked" exceptions and rethrow them as 
	        // runtime exceptions, which are the only kind of exceptions static 
	        // blocks are allowed to throw. See 
            // <http://stackoverflow.com/a/15289277/402891>
	        throw new ExceptionInInitializerError(e);
	    }
        
    }
}
