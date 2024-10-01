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
}