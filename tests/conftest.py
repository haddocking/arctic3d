import pytest


@pytest.fixture
def compare_responses():
    def are_keys_present(sub_dict: dict, main_dict: dict) -> bool:
        """
        Check if all keys in 'sub_dict' are present in 'main_dict', including nested dictionaries.
        This function is recursive and handles any depth of nesting.
        """
        if not isinstance(sub_dict, dict) or not isinstance(main_dict, dict):
            return False

        for key, value in sub_dict.items():
            if key in main_dict:
                if isinstance(value, dict):
                    # Recursively check nested dictionaries.
                    if not isinstance(
                        main_dict[key], dict
                    ) or not are_keys_present(value, main_dict[key]):
                        return False
                # If the value is not a dict, we only care about the key which we know is present
            else:
                # If a key is not present at any level, return False
                return False
        return True

    return are_keys_present
