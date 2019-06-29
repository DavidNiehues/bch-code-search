import math as m
import logging

logging.basicConfig()
logger = logging.getLogger("find_good_bch_codes")
logger.setLevel(logging.WARNING)


def list_suitable_bch_codes(n, min_dimension, required_relative_minimal_bose_distance):
    # the minimal distance required in order to achieve the given relative distance.
    required_min_distance = m.ceil(n * required_relative_minimal_bose_distance)

    # compute all 2-cyclotomic cosets mod n.
    cosets = q_cyclotomic_cosets_mod_n(2, n)
    print("Finished computing cyclotomic cosets")

    smallest_n = n + 1
    best_k = 0
    best_d = 0
    best_puncturing = 0

    # consider all possible starting points for sequences for 2-cyclotomic cosets mod n
    for start in range(0, n - 1):
        # initialize the defining set to be empty and i to be at the start position.
        defining_set = set()
        i = start

        # add 2-cyclotomic cosets mod n to the defining set until the required minimal distance is reached or there are
        # no further 2-cyclotomic cosets mod n.
        while i <= n - 1:
            defining_set = defining_set.union(cosets[i])
            i += 1
            achieved_bose_distance = _number_of_consecutive_elements(defining_set)

            # the required minimal distance is achieved if there are more than required_min_distance consecutive elements in
            # the defining set.
            if achieved_bose_distance >= required_min_distance:

                # the dimension of the code is n minus the size of the defining set (Theorem 4.4.2 Huffmann-Pless)
                dimension = n - len(defining_set)
                if dimension > min_dimension:
                    num_puncturing = _get_puncturing_amount(achieved_bose_distance,
                                                            required_relative_minimal_bose_distance,
                                                            n)

                    logger.debug("[{:d}, {:d}, {:d}]_2 at start {:d} and for {:d} cosets".format(n, dimension,
                                                                                                 achieved_bose_distance,
                                                                                                 start,
                                                                                                 i - start))
                    if n - num_puncturing < smallest_n:
                        smallest_n = n - num_puncturing
                        best_puncturing = num_puncturing
                        best_k = dimension
                        best_d = achieved_bose_distance

                        print("[{:d}, {:d}, {:d}] => {:d} x puncturing => [{:d}, {:d}, {:d}]".format(n, dimension,
                                                                                                     achieved_bose_distance,
                                                                                                     num_puncturing,
                                                                                                     n - num_puncturing,
                                                                                                     dimension,
                                                                                                     achieved_bose_distance
                                                                                                     - num_puncturing))

    print(
        "Best: [{:d}, {:d}, {:d}] => {:d} x puncturing => [{:d}, {:d}, {:d}]".format(n, best_k, best_d, best_puncturing,
                                                                                     n - best_puncturing, best_k,
                                                                                     best_d - best_puncturing))


def q_cyclotomic_cosets_mod_n(q, n):
    """ Computes all q-cyclotomic cosets modulo n.
		The result is a dictionary mapping all 0 <= i <= n
		to the cyclotomic coset containing i
	"""
    cosets = dict()
    zero_coset = set()
    zero_coset.add(0)
    cosets[0] = zero_coset

    for i in range(1, n):
        if not i in cosets:
            coset = set()
            coset.add(i)
            cosets[i] = coset
            next_element = i * q % n
            while not next_element == i:
                coset.add(next_element)
                cosets[next_element] = coset
                next_element = next_element * q % n
    return cosets


def _number_of_consecutive_elements(a):
    """ Given the set a, the function outputs the largest number of consecutive elements in a
    """

    # A set of size zero always has zero consectutive numbers in it.
    if len(a) == 0:
        return 0

    sorted_elements = sorted(list(a))

    # the length of the longest sequence found so far
    longest_sequence_length = 1
    # the index where the longest sequence found so far starts (in the sorted list)
    longest_sequence_start = 0
    # the length of the currently considered sequence
    current_sequence_length = 1
    # the index where the currently considered sequence starts
    current_sequence_start = 0

    # initialize the previous element with the first element in order to emulate a do-while loop
    previous_element = sorted_elements[0]
    for i in range(1, len(sorted_elements)):
        current_element = sorted_elements[i]
        if current_element == previous_element + 1:
            # if the current sequence is continued, update the length of the current
            # sequence and, if it is longer than the currently longest sequence, update the information
            # on the longest sequence found so far to the current sequence.
            current_sequence_length += 1
            if current_sequence_length > longest_sequence_length:
                longest_sequence_length = current_sequence_length
                longest_sequence_start = current_sequence_start
        else:
            # if the current element is not the continuation of
            # the previous sequence, start a new sequence of length 1
            current_sequence_start = i
            current_sequence_length = 1
        previous_element = current_element
    return longest_sequence_length


def _get_puncturing_amount(achieved_bose_distance, required_rel_min_distance, code_out_length):
    """ Given the achieved bose distance, the output length of an error correcting code and a desired relative minimal
     distance, this function computes how often the code can be punctured while still maintaining the desired relative
      minimal distance.
    """
    return m.floor(
        (achieved_bose_distance - code_out_length * required_rel_min_distance) / (1 - required_rel_min_distance))
