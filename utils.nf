// utils.nf

def checkMeta(item, expectedMeta = []) {
    // String 
    // NullObject
    // String[]: List of Strings
    // ArrayList: Empty list
    def meta = item[0]
    if (!meta.containsKey("id")) {
        throw new IllegalArgumentException("Missing metadata key: id")
    }

    expectedMeta.each { key, expectedType ->
        
        if (!meta.containsKey(key)) {
            throw new IllegalArgumentException("Missing metadata key: ${key}")
        }
        
        def value = meta.get(key)

        def actualType = value.getClass().getSimpleName()

        if (!expectedType.contains(actualType)) {
            throw new IllegalArgumentException("For item '${meta.id}', metadata key '${key}' is of type '${actualType}', expected '${expectedType}'")
        }
    }

    def checkFiles
    checkFiles = { value ->
        if (value == null) {
        } else if (value instanceof Path || value instanceof File) {
            if (!value.exists()) {
                throw new FileNotFoundException("File '${value}' does not exist")
            }
        } else if (value instanceof Collection || value.getClass().isArray()) {
            value.each { v ->
                checkFiles(v)
            }
        }
    }

    item[1..-1].each { val ->
        checkFiles(val)
    }
}

