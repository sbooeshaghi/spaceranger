# Based on the LZ-String javascript found here:
#   http://pieroxy.net/blog/pages/lz-string/index.html
#   version 1.4.4

keyStrUriSafe = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-$"


def compressToEncodedURIComponent(input_data):
    if input_data is None:
        return ""

    return _compress(input_data, 6, lambda a: keyStrUriSafe[a])


def _compress(uncompressed, bits_per_char, get_char_from_int):
    #pylint: disable=too-many-branches,too-many-statements
    if uncompressed is None:
        return ""

    context_dictionary = {}
    context_dictionaryToCreate = {}
    context_wc = ""
    context_w = ""
    context_enlargeIn = 2  # Compensate for the first entry which should not count
    context_dictSize = 3
    context_numBits = 2
    context_data = []
    context_data_val = 0
    context_data_position = 0

    def get_char_iter(s, chunk_size=1024):
        if hasattr(s, '__len__'):
            for c in s:
                yield c
        elif hasattr(s, 'read'):
            while True:
                chunk = s.read(chunk_size)
                if len(chunk) == 0:
                    break
                for c in chunk:
                    yield c
        else:
            raise Exception(
                "Don't know how to compress object of type %s" % type(s))

    for context_c in get_char_iter(uncompressed):
        if context_c not in context_dictionary:
            context_dictionary[context_c] = context_dictSize
            context_dictSize += 1
            context_dictionaryToCreate[context_c] = True

        context_wc = context_w + context_c
        if context_wc in context_dictionary:
            context_w = context_wc
        else:
            if context_w in context_dictionaryToCreate:
                if ord(context_w[0]) < 256:
                    for _ in xrange(context_numBits):
                        context_data_val <<= 1

                        if context_data_position == bits_per_char - 1:
                            context_data_position = 0
                            context_data.append(
                                get_char_from_int(context_data_val))
                            context_data_val = 0
                        else:
                            context_data_position += 1

                    value = ord(context_w[0])
                    for _ in xrange(8):
                        context_data_val = (
                            context_data_val << 1) | (value & 1)
                        if context_data_position == bits_per_char - 1:
                            context_data_position = 0
                            context_data.append(
                                get_char_from_int(context_data_val))
                            context_data_val = 0
                        else:
                            context_data_position += 1

                        value >>= 1
                else:
                    value = 1
                    for _ in xrange(context_numBits):
                        context_data_val = (context_data_val << 1) | value
                        if context_data_position == bits_per_char - 1:
                            context_data_position = 0
                            context_data.append(
                                get_char_from_int(context_data_val))
                            context_data_val = 0
                        else:
                            context_data_position += 1

                        value = 0

                    value = ord(context_w[0])
                    for _ in xrange(16):
                        context_data_val = (
                            context_data_val << 1) | (value & 1)
                        if context_data_position == bits_per_char - 1:
                            context_data_position = 0
                            context_data.append(
                                get_char_from_int(context_data_val))
                            context_data_val = 0
                        else:
                            context_data_position += 1

                        value >>= 1

                context_enlargeIn -= 1
                if context_enlargeIn == 0:
                    context_enlargeIn = 2**context_numBits
                    context_numBits += 1
                context_dictionaryToCreate.pop(context_w, None)
            else:
                value = context_dictionary[context_w]

                for _ in xrange(context_numBits):
                    context_data_val = (context_data_val << 1) | (value & 1)
                    if context_data_position == bits_per_char - 1:
                        context_data_position = 0
                        context_data.append(
                            get_char_from_int(context_data_val))
                        context_data_val = 0
                    else:
                        context_data_position += 1

                    value >>= 1

            context_enlargeIn -= 1
            if context_enlargeIn == 0:
                context_enlargeIn = 2**context_numBits
                context_numBits += 1

            # Add wc to the dictionary.
            context_dictionary[context_wc] = context_dictSize
            context_dictSize += 1
            context_w = context_c

    # Output the code for w.
    if context_w != '':
        if context_w in context_dictionaryToCreate:
            if ord(context_w[0]) < 256:
                for _ in xrange(context_numBits):
                    context_data_val <<= 1
                    if context_data_position == bits_per_char - 1:
                        context_data_position = 0
                        context_data.append(
                            get_char_from_int(context_data_val))
                        context_data_val = 0
                    else:
                        context_data_position += 1

                value = ord(context_w[0])
                for _ in xrange(8):
                    context_data_val = (context_data_val << 1) | (value & 1)
                    if context_data_position == bits_per_char - 1:
                        context_data_position = 0
                        context_data.append(
                            get_char_from_int(context_data_val))
                        context_data_val = 0
                    else:
                        context_data_position += 1

                    value >>= 1
            else:
                value = 1
                for _ in xrange(context_numBits):
                    context_data_val = (context_data_val << 1) | value
                    if context_data_position == bits_per_char - 1:
                        context_data_position = 0
                        context_data.append(
                            get_char_from_int(context_data_val))
                        context_data_val = 0
                    else:
                        context_data_position += 1

                    value = 0

                value = ord(context_w[0])
                for _ in xrange(16):
                    context_data_val = (context_data_val << 1) | (value & 1)
                    if context_data_position == bits_per_char - 1:
                        context_data_position = 0
                        context_data.append(
                            get_char_from_int(context_data_val))
                        context_data_val = 0
                    else:
                        context_data_position += 1

                    value >>= 1

            context_enlargeIn -= 1
            if context_enlargeIn == 0:
                context_enlargeIn = 2**context_numBits
                context_numBits += 1
            context_dictionaryToCreate.pop(context_w, None)
        else:
            value = context_dictionary[context_w]

            for _ in xrange(context_numBits):
                context_data_val = (context_data_val << 1) | (value & 1)
                if context_data_position == bits_per_char - 1:
                    context_data_position = 0
                    context_data.append(get_char_from_int(context_data_val))
                    context_data_val = 0
                else:
                    context_data_position += 1

                value >>= 1

        context_enlargeIn -= 1
        if context_enlargeIn == 0:
            context_enlargeIn = 2**context_numBits
            context_numBits += 1

    # Mark the end of the stream
    value = 2
    for _ in xrange(context_numBits):
        context_data_val = (context_data_val << 1) | (value & 1)
        if context_data_position == bits_per_char - 1:
            context_data_position = 0
            context_data.append(get_char_from_int(context_data_val))
            context_data_val = 0
        else:
            context_data_position += 1

        value >>= 1

    # Flush the last char
    while True:
        context_data_val <<= 1
        if context_data_position == bits_per_char - 1:
            context_data.append(get_char_from_int(context_data_val))
            break
        else:
            context_data_position += 1

    return ''.join(context_data)
